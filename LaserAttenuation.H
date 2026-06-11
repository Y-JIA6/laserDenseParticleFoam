/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  10
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

// =============================================================================
// LaserAttenuation.H  —  TEM00 Gaussian Laser Beam Intensity Field
// =============================================================================
//
// PURPOSE
// -------
//   Compute the Eulerian volumetric field  laserIntensity [W/m²]  at every
//   mesh cell centre for a TEM00 (fundamental-mode) Gaussian beam.
//   The field is used for visualisation; particle energy absorption is
//   computed independently (analytically) in particleInfo.H.
//
// PHYSICS: TEM00 Gaussian Beam
// ----------------------------
//   The TEM00 mode has a Gaussian transverse intensity profile.  In cylindrical
//   coordinates (x along the beam axis, r radial):
//
//     I(x, r) = (2P / (pi * Rb(x)^2)) * exp(-2*r^2 / Rb(x)^2)       ...(1)
//
//   This form satisfies the power normalisation exactly:
//
//     integral_0^inf  I(x,r) * 2*pi*r dr
//     = (2P / Rb^2) * integral_0^inf exp(-2*r^2/Rb^2) * r dr
//     = (2P / Rb^2) * [Rb^2 / 4]  *  4  *  (1/2)
//     = P                                                              ...(2)
//
//   i.e. the total power passing through any cross-section equals P [W],
//   independent of Rb.  The factor 2 (not the common 1/e² factor of 1) comes
//   from expressing the exponent as -2r²/Rb², so Rb is the 1/e² radius.
//
//   Peak intensity on axis (r=0):
//
//     I_peak(x) = 2P / (pi * Rb(x)^2)                                 ...(3)
//
// BEAM RADIUS MODEL  Rb(x)
// -------------------------
//   Real focussed laser heads have three zones along x:
//
//   Zone 1 — converging  (x < xf - DOF/2):
//     The beam converges toward the focal point from the upstream side.
//
//   Zone 2 — flat focus  (xf - DOF/2  <=  x  <=  xf + DOF/2):
//     Within the Depth-of-Field (DOF) the beam radius is approximately
//     constant at the minimum value R0 (diffraction-limited waist).
//
//   Zone 3 — diverging   (x > xf + DOF/2):
//     The beam diverges at half-angle theta after leaving the focal zone.
//
//   Compact unified formula:
//
//     dist_to_DOF(x) = max(0, xf - DOF/2 - x)     [upstream gap]
//                    + max(0, x - xf - DOF/2)      [downstream gap]
//
//     Rb(x) = R0 + dist_to_DOF(x) * tan(theta)                        ...(4)
//
//   Parameters from LaserProperties dictionary:
//     laserPower      P       [W]    Total beam power
//     beamRadiusMin   R0      [m]    Minimum (focal) 1/e² beam radius
//     halfAngle       theta   [rad]  Convergence/divergence half-angle
//     focalDistance   xf      [m]    Axial position of focal point
//     depthOfField    DOF     [m]    Length of flat-focus zone
//
// MOVING LASER HEAD
// -----------------
//   The beam origin advances with time:
//
//     laserOrigin(t) = laserInitialPosition + laserVelocity * t        ...(5)
//
//   When  laserVelocity == 0  the field is time-invariant and is computed
//   once only (static cache flag laserIntensityComputed).
//
// TRUNCATION OPTIMISATION
// -----------------------
//   For  r > 3*Rb  the Gaussian factor  exp(-2*(3)^2) = exp(-18) ≈ 1.5e-8,
//   contributing < 2e-6 % of peak intensity.  These cells are set to zero,
//   saving ~97% of exp() evaluations in typical meshes.
//
// OUTPUT
// ------
//   laserIntensity  [W/m²]   volScalarField, written at every outputTime.
// =============================================================================

{
    // =========================================================================
    // Read LaserProperties
    // =========================================================================
    const IOdictionary LaserProperties
    (
        IOobject
        (
            "LaserProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar pi = Foam::constant::mathematical::pi;

    // ---- Beam direction (normalised) ----------------------------------------
    vector laserDir = LaserProperties.get<vector>("laserDirection");
    laserDir /= mag(laserDir);

    // ---- Laser origin (supports moving head: origin = init + vel*t) ---------
    vector laserOriginBase = LaserProperties.get<vector>("laserOrigin");
    vector laserInitPos    = LaserProperties.getOrDefault<vector>
                             ("laserInitialPosition", laserOriginBase);
    vector laserVel        = LaserProperties.getOrDefault<vector>
                             ("laserVelocity", vector::zero);
    const vector laserOrigin = laserInitPos + laserVel * runTime.value();

    // ---- Beam profile parameters --------------------------------------------
    const scalar P     = LaserProperties.get<scalar>("laserPower");     // [W]
    const scalar R0    = LaserProperties.get<scalar>("beamRadiusMin");  // [m]
    const scalar theta = LaserProperties.get<scalar>("halfAngle");      // [rad]
    const scalar xf    = LaserProperties.get<scalar>("focalDistance");  // [m]
    const scalar DOF   = LaserProperties.get<scalar>("depthOfField");   // [m]

    const scalar xDOF_lo  = xf - 0.5*DOF;   // upstream boundary of flat-focus zone
    const scalar xDOF_hi  = xf + 0.5*DOF;   // downstream boundary of flat-focus zone
    const scalar tanTheta = Foam::tan(theta);

    // =========================================================================
    // Static laser optimisation
    // If the laser is stationary, I(x,r) is time-invariant: compute once only.
    // =========================================================================
    const bool laserIsMoving = (mag(laserVel) > SMALL);

    // static bool persists across all time-step calls of this #include block
    static bool laserIntensityComputed = false;

    if (laserIsMoving || !laserIntensityComputed)
    {
        Info<< "LaserAttenuation: computing intensity field at t = "
            << runTime.timeName() << nl
            << "  laserDirection = " << laserDir     << nl
            << "  laserOrigin    = " << laserOrigin  << nl
            << "  laserPower     = " << P            << " W" << nl
            << "  beamRadiusMin  = " << R0*1e3       << " mm" << nl
            << "  halfAngle      = " << theta        << " rad" << nl
            << "  focalDistance  = " << xf*1e3       << " mm" << nl
            << "  depthOfField   = " << DOF*1e3      << " mm" << endl;

        // ---------------------------------------------------------------------
        // Eq.(4): distance from cell centre to the nearest DOF boundary,
        //         zero when the cell is inside the flat-focus zone.
        // Eq.(4): Rb = R0 + dist_to_DOF * tan(theta)
        // ---------------------------------------------------------------------
        forAll(mesh.C(), celli)
        {
            const vector relPos = mesh.C()[celli] - laserOrigin;

            const scalar x = relPos & laserDir;                 // axial distance

            // x < 0: cell is upstream of laserOrigin (inside nozzle focusing optics).
            // The free-space Gaussian model is only valid downstream of the nozzle
            // exit.  Setting I = 0 here prevents a spurious low-intensity tail from
            // extending into the nozzle gas channels and matches the paper's domain.
            if (x < 0)
            {
                laserIntensity[celli] = 0.0;
                continue;
            }

            const scalar distToDOF = max(scalar(0), xDOF_lo - x)
                                   + max(scalar(0), x - xDOF_hi);
            const scalar Rb = R0 + distToDOF*tanTheta;

            const scalar r = mag(relPos - x*laserDir);         // radial distance

            // Truncation: exp(-18) < 1.5e-8 — negligible contribution
            if (r > 3.0*Rb)
            {
                laserIntensity[celli] = 0.0;
            }
            else
            {
                // Eq.(1): I = 2P/(pi*Rb^2) * exp(-2*r^2/Rb^2)
                laserIntensity[celli] =
                    (2.0*P / (pi*sqr(Rb))) * Foam::exp(-2.0*sqr(r)/sqr(Rb));
            }
        }
        laserIntensity.correctBoundaryConditions();

        Info<< "  Peak laser intensity = "
            << gMax(laserIntensity) << " W/m^2" << endl;

        laserIntensityComputed = true;
    }
    else
    {
        Info<< "LaserAttenuation: static laser - reusing pre-computed field"
            << " (t = " << runTime.timeName() << ")" << endl;
    }

    // =========================================================================
    // Write laserIntensity at every output time
    // =========================================================================
    if (runTime.outputTime())
    {
        laserIntensity.write();
    }
}

