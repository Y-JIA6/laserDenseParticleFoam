/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  2506
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

// =============================================================================
// particleInfo.H  —  Per-Particle Thermal Energy Balance
// =============================================================================
//
// PURPOSE
// -------
//   Advance each particle's temperature by one CFD time step dt using an
//   explicit Euler integration of the enthalpy ODE:
//
//     m * dH/dt = Q_laser - Q_rad - Q_conv                             ...(A)
//
//   where H is the specific enthalpy [J/kg] defined to include latent heat.
//
// -----------------------------------------------------------------------------
// 1. LASER ENERGY INPUT  Q_laser  [W]
// -----------------------------------------------------------------------------
//
//   The particle intercepts the laser beam over its projected cross-section:
//
//     Q_laser = alpha_laser * I(x_p, r_p) * A_proj                    ...(1)
//
//   where:
//     alpha_laser  — laser absorptivity of the particle material [-]
//                    (wavelength-dependent; 0.47 for Stellite 6 at ~1 µm)
//     I(x_p, r_p) — TEM00 Gaussian intensity [W/m²] at the particle location,
//                    computed analytically (see below; same formula as
//                    LaserAttenuation.H, equation 1)
//     A_proj       = pi/4 * d^2   particle projected area [m²]
//
//   Gaussian intensity at particle position (x_p along beam axis, r_p radial):
//
//     I = (2P / (pi * Rb^2)) * exp(-2 * r_p^2 / Rb^2)                 ...(2)
//     Rb = R0 + dist_to_DOF(x_p) * tan(theta)                          ...(3)
//
//   Note: I is evaluated analytically per particle, NOT by interpolating
//   the laserIntensity Eulerian field, to avoid mesh-resolution errors.
//
//   Truncation: when r_p > 3*Rb, exp(-18) < 1.5e-8 → Q_laser set to 0.
//
// -----------------------------------------------------------------------------
// 2. RADIATIVE LOSS  Q_rad  [W]
// -----------------------------------------------------------------------------
//
//   Stefan-Boltzmann radiation from the full particle surface to the ambient:
//
//     Q_rad = epsilon * sigma_SB * A_surf * (T_p^4 - T_amb^4)         ...(4)
//
//   where:
//     epsilon    — hemispherical total emissivity [-]
//     sigma_SB   = 5.67e-8  W/(m²·K⁴)   Stefan-Boltzmann constant
//     A_surf     = pi * d^2              full sphere surface area [m²]
//     T_p        — current particle temperature [K]
//     T_amb      — ambient (surroundings) temperature [K]
//
//   Q_rad > 0 when T_p > T_amb (heat loss); Q_rad < 0 when T_p < T_amb
//   (radiative gain, rare in practice).
//
// -----------------------------------------------------------------------------
// 3. CONVECTIVE LOSS  Q_conv  [W]
// -----------------------------------------------------------------------------
//
//   Forced convection using the Ranz-Marshall (1952) correlation:
//
//     Q_conv = h * A_surf * (T_p - T_amb)                              ...(5)
//     h      = Nu * k_g / d                                            ...(6)
//     Nu     = 2 + 0.6 * Re_p^(1/2) * Pr_g^(1/3)      [Ranz-Marshall] ...(7)
//     Re_p   = rho_c * |U_c - U_p| * d / mu_c                         ...(8)
//
//   where:
//     k_g    — gas thermal conductivity [W/(m·K)]
//     Pr_g   — gas Prandtl number [-]  (Ar: 0.741)
//     Re_p   — particle Reynolds number [-]
//     rho_c, mu_c, U_c  — local gas density, viscosity, velocity
//                          queried from the CFD cell containing the particle
//
//   Limiting case: if p.cell() < 0 (particle outside mesh), Nu = 2 → h = 2*k_g/d
//   (pure conduction in a quiescent medium, Stokes limit).
//
// -----------------------------------------------------------------------------
// 4. ENTHALPY METHOD FOR PHASE CHANGE
// -----------------------------------------------------------------------------
//
//   Direct integration of  dT/dt = (Q_net)/(m*cp)  fails across the mushy
//   zone because cp is formally infinite at melting.  The enthalpy method
//   avoids this by integrating H instead:
//
//     H(T) = cp * T                                       (solid,  T <= T_sol)
//     H(T) = cp * T + L_f*(T - T_sol)/(T_liq - T_sol)   (mushy zone)      ...(9)
//     H(T) = cp * T + L_f                                (liquid, T >= T_liq)
//
//   Explicit Euler step:
//
//     H_new = H_old + (Q_laser - Q_rad - Q_conv) / m * dt              ...(10)
//
//   Inverse mapping H → T:
//
//     T = H_new / cp                           if H_new <= H_sol = cp*T_sol
//     T = (H_new*dTL + L_f*T_sol)/(cp*dTL+Lf) if H_sol < H_new < H_liq  ...(11)
//     T = (H_new - L_f) / cp                   if H_new >= H_liq = cp*T_liq+Lf
//
//   where  dTL = T_liq - T_sol.
//
//   The mushy-zone formula in (11) is derived by solving (9) for T:
//     H = cp*T + L_f*(T - T_sol)/dTL
//     H*dTL = cp*T*dTL + L_f*T - L_f*T_sol
//     T*(cp*dTL + L_f) = H*dTL + L_f*T_sol
//     T = (H*dTL + L_f*T_sol) / (cp*dTL + L_f)
//
// -----------------------------------------------------------------------------
// 5. EXPLICIT EULER STABILITY CONDITION
// -----------------------------------------------------------------------------
//
//   For explicit Euler to remain stable the thermal time constant must satisfy:
//
//     dt  <=  2 * tau_thermal                                          ...(12)
//     tau_thermal = m*cp / beta_eff
//     beta_eff    = A_surf * (h_conv + 4*epsilon*sigma_SB*T_p^3)      ...(13)
//
//   The 4T^3 term is the linearised radiation coefficient (first-order Taylor
//   expansion of T^4 about T_p).  If (12) is violated a WarningInFunction is
//   issued; consider reducing deltaT or using an implicit ODE solver.
//
// MPI / PARALLELISM
// -----------------
//   Temperature is stored in a static  Map<scalar> T_persistent  keyed by
//   (origProc * 1e6 + origId).  After each step all procs exchange their
//   freshly-computed values via gatherList/broadcastList so that the current
//   owner's value always supersedes any stale copies left from particle
//   migration between processors.
//
// =============================================================================

{
    // -----------------------------------------------------------------------
    // Read material and laser properties (statics, read once per run)
    // -----------------------------------------------------------------------
    static bool propsRead = false;
    // Shared (gas/environment) properties — same for all materials
    static scalar sigma_SB, T_amb, k_g, Pr_g, Pr_g_13;
    // Per-material properties — keyed by injectorID (= p.typeId())
    // injectorID -1 (unset) falls back to key 0
    static Map<scalar> mat_alpha, mat_epsilon, mat_cp, mat_rhoP;
    static Map<scalar> mat_Lf, mat_Tsol, mat_Tliq, mat_Tinit;
    static scalar laserP, laserR0, laserTheta, laserXf, laserDOF;
    static vector laserInitPos, laserVel, laserDir;

    if (!propsRead)
    {
        const IOdictionary particleProperties
        (
            IOobject
            (
                "particleProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // --- Shared environment/gas properties ---
        sigma_SB = particleProperties.getOrDefault<scalar>("sigma", 5.67e-8);
        T_amb    = particleProperties.getOrDefault<scalar>("T_amb", 300.0);
        k_g      = particleProperties.getOrDefault<scalar>("k_g",   0.01580);
        Pr_g     = particleProperties.getOrDefault<scalar>("Pr_g",  0.741);
        Pr_g_13  = Foam::pow(Pr_g, 1.0/3.0);

        // --- Per-material properties ---
        // Multi-material mode: read 'materials' list (index = injectorID/typeId)
        // Single-material mode: read flat scalars, wrap into index-0 (backward compatible)
        if (particleProperties.found("materials"))
        {
            const dictionary& matsDict = particleProperties.subDict("materials");
            const wordList matNames(matsDict.toc());

            Info<< "particleInfo: multi-material mode (" << matNames.size() << " materials)" << nl
                << "  sigma_SB = " << sigma_SB << " W/(m^2·K^4)" << nl
                << "  T_amb    = " << T_amb    << " K"            << nl
                << "  k_g      = " << k_g      << " W/(m·K)"      << nl
                << "  Pr_g     = " << Pr_g     << " [-]"          << endl;

            forAll(matNames, i)
            {
                const dictionary& d = matsDict.subDict(matNames[i]);
                // injectorID in the material dict is the key for the Map
                const label id = d.get<label>("injectorID");
                mat_alpha.set  (id, d.get<scalar>("alpha_laser"));
                mat_epsilon.set(id, d.get<scalar>("epsilon"));
                mat_cp.set     (id, d.get<scalar>("cp"));
                mat_rhoP.set   (id, d.get<scalar>("rhoP"));
                mat_Lf.set     (id, d.get<scalar>("L_f"));
                mat_Tsol.set   (id, d.get<scalar>("T_sol"));
                mat_Tliq.set   (id, d.get<scalar>("T_liq"));
                mat_Tinit.set  (id, d.getOrDefault<scalar>("T_init", T_amb));
                Info<< "  Material '" << matNames[i] << "' (injectorID=" << id << "):" << nl
                    << "    rhoP=" << mat_rhoP[id]  << " kg/m^3"  << nl
                    << "    cp="   << mat_cp[id]    << " J/(kg·K)" << nl
                    << "    alpha_laser=" << mat_alpha[id]   << " [-]" << nl
                    << "    epsilon="     << mat_epsilon[id] << " [-]" << nl
                    << "    L_f="  << mat_Lf[id]   << " J/kg" << nl
                    << "    T_sol="<< mat_Tsol[id] << " K  T_liq=" << mat_Tliq[id] << " K" << nl
                    << "    T_init=" << mat_Tinit[id] << " K" << endl;
            }
        }
        else
        {
            // Single-material (backward compatible): store at key 0
            mat_alpha.set  (0, particleProperties.getOrDefault<scalar>("alpha_laser", 0.4));
            mat_epsilon.set(0, particleProperties.getOrDefault<scalar>("epsilon",     0.3));
            mat_cp.set     (0, particleProperties.getOrDefault<scalar>("cp",          421.0));
            mat_rhoP.set   (0, particleProperties.getOrDefault<scalar>("rhoP",        8380.0));
            mat_Lf.set     (0, particleProperties.getOrDefault<scalar>("L_f",         2.27e5));
            mat_Tsol.set   (0, particleProperties.getOrDefault<scalar>("T_sol",       1558.0));
            mat_Tliq.set   (0, particleProperties.getOrDefault<scalar>("T_liq",       1668.0));
            mat_Tinit.set  (0, particleProperties.getOrDefault<scalar>("T_init",      T_amb));

            Info<< "particleInfo: single-material mode" << nl
                << "  rhoP        = " << mat_rhoP[0]    << " kg/m^3"       << nl
                << "  cp          = " << mat_cp[0]      << " J/(kg·K)"     << nl
                << "  alpha_laser = " << mat_alpha[0]   << " [-]"          << nl
                << "  epsilon     = " << mat_epsilon[0] << " [-]"          << nl
                << "  L_f         = " << mat_Lf[0]      << " J/kg"         << nl
                << "  T_sol       = " << mat_Tsol[0]    << " K"            << nl
                << "  T_liq       = " << mat_Tliq[0]    << " K"            << nl
                << "  T_amb       = " << T_amb           << " K"            << nl
                << "  T_init      = " << mat_Tinit[0]   << " K"            << nl
                << "  k_g         = " << k_g             << " W/(m·K)"     << nl
                << "  Pr_g        = " << Pr_g            << " [-]"         << nl
                << "  sigma_SB    = " << sigma_SB        << " W/(m^2·K^4)" << endl;
        }
        propsRead = true;
    }

    // Read laser geometry every step (supports moving laser head)
    {
        const IOdictionary laserProps
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

        laserP     = laserProps.get<scalar>("laserPower");
        laserR0    = laserProps.get<scalar>("beamRadiusMin");
        laserTheta = laserProps.get<scalar>("halfAngle");
        laserXf    = laserProps.get<scalar>("focalDistance");
        laserDOF   = laserProps.get<scalar>("depthOfField");

        vector laserOriginBase = laserProps.get<vector>("laserOrigin");
        laserInitPos = laserProps.getOrDefault<vector>("laserInitialPosition", laserOriginBase);
        laserVel     = laserProps.getOrDefault<vector>("laserVelocity", vector::zero);
        laserDir     = laserProps.get<vector>("laserDirection");
        laserDir    /= mag(laserDir);
    }

    const scalar pi = Foam::constant::mathematical::pi;
    const scalar dt = runTime.deltaTValue();

    // -----------------------------------------------------------------------
    // Persistent temperature map -- key = origProc*1e6 + origId
    // Cleared only at solver restarts (static lifetime).
    // -----------------------------------------------------------------------
    static Map<scalar> T_persistent;

    // -----------------------------------------------------------------------
    // Compute laser beam profile at current time
    // -----------------------------------------------------------------------
    const vector laserOriginNow = laserInitPos + laserVel * runTime.value();
    const scalar xDOF_lo  = laserXf - 0.5*laserDOF;
    const scalar xDOF_hi  = laserXf + 0.5*laserDOF;
    const scalar tanTheta = Foam::tan(laserTheta);

    // -----------------------------------------------------------------------
    // Step 1: Compute new T for every locally-owned (active) particle
    // -----------------------------------------------------------------------
    label nOwned = kinematicCloud.size();
    labelList ownedKeys(nOwned);
    scalarList ownedTemps(nOwned);
    label owned = 0;

    label nHeated = 0;
    scalar T_max_local   = -GREAT;
    scalar dT_max_local  =  0.0;
    vector pos_dTmax_local = vector::zero;
    vector pos_Tmax_local  = vector::zero;

    forAllConstIter(basicKinematicTypeCloud, kinematicCloud, iter)
    {
        const auto& p = iter();
        const label key = p.origProc() * 1000000 + p.origId();

        // Material lookup by injectorID (= p.typeId()).
        // typeId=-1 means injectorID not set in kinematicCloudProperties → fall back to key 0.
        const label matId = (p.typeId() >= 0 && mat_alpha.found(p.typeId()))
                          ? p.typeId() : 0;
        const scalar alpha_laser = mat_alpha[matId];
        const scalar epsilon_p   = mat_epsilon[matId];
        const scalar cp_p        = mat_cp[matId];
        const scalar rhoP        = mat_rhoP[matId];
        const scalar L_f         = mat_Lf[matId];
        const scalar T_sol       = mat_Tsol[matId];
        const scalar T_liq       = mat_Tliq[matId];
        const scalar T_init_mat  = mat_Tinit[matId];

        // Retrieve temperature from map; T_init_mat for brand-new particles
        const scalar T_old = T_persistent.lookup(key, T_init_mat);
        scalar T_curr = T_old;

        // --- Particle geometry ---
        const scalar diam   = p.d();
        const scalar A_proj = pi / 4.0 * sqr(diam);
        const scalar A_surf = pi * sqr(diam);
        const scalar V      = pi / 6.0 * pow3(diam);
        const scalar m      = rhoP * V;

        // --- Laser intensity: analytic TEM00 Gaussian (Eqs. 2-3) ---
        // x_l: axial projection of (particle - laserOrigin) onto beam axis
        // r_l: radial distance from beam axis
        // Rb:  effective beam radius at this axial position
        const vector relPos  = p.position() - laserOriginNow;
        const scalar x_l     = relPos & laserDir;
        const scalar r_l     = mag(relPos - x_l * laserDir);

        scalar I_laser = 0.0;
        // x_l < 0: particle is upstream of laserOrigin (inside nozzle optics).
        // The free-space Gaussian is only valid downstream of the nozzle exit.
        if (x_l >= 0)
        {
            const scalar distDOF = max(scalar(0), xDOF_lo - x_l)
                                 + max(scalar(0), x_l - xDOF_hi);
            const scalar Rb      = laserR0 + distDOF * tanTheta;

            if (r_l <= 3.0 * Rb)       // truncate below ~1.5e-8 of peak
            {
                I_laser = (2.0 * laserP / (pi * sqr(Rb)))
                        * Foam::exp(-2.0 * sqr(r_l) / sqr(Rb));
            }
        }

        // Eq.(1): Q_laser = alpha * I * A_proj
        const scalar Q_laser = alpha_laser * I_laser * A_proj;

        // Eq.(4): Stefan-Boltzmann radiation  Q_rad = eps*sigma*(T^4 - T_amb^4)
        const scalar Q_rad = epsilon_p * sigma_SB * A_surf
                           * (pow4(T_curr) - pow4(T_amb));

        // Eqs.(5-8): Ranz-Marshall convection  Q_conv = h*A_surf*(T_p - T_amb)
        const label cellI = p.cell();
        scalar h_conv = 2.0 * k_g / diam;   // Stokes limit Nu=2 (fallback)
        if (cellI >= 0)
        {
            const scalar rhoLoc = rhoc[cellI];
            const scalar muLoc  = muc[cellI];
            const scalar Re_p   = rhoLoc * mag(Uc[cellI] - p.U()) * diam
                                / max(muLoc, SMALL);
            const scalar Nu     = 2.0 + 0.6 * Foam::sqrt(Re_p) * Pr_g_13;
            h_conv = Nu * k_g / diam;
        }
        const scalar Q_conv = h_conv * A_surf * (T_curr - T_amb);

        // Eqs.(9-11): enthalpy method — advance H then invert to T
        //
        // Boundary enthalpies:
        //   H_sol = cp * T_sol                (top of solid range)
        //   H_liq = cp * T_liq + L_f          (bottom of liquid range)
        //
        // H → T inversion:
        //   solid  : T = H / cp
        //   mushy  : T = (H*dTL + L_f*T_sol) / (cp*dTL + L_f)   [Eq.(11)]
        //   liquid : T = (H - L_f) / cp
        const scalar H_sol = cp_p * T_sol;
        const scalar H_liq = cp_p * T_liq + L_f;
        const scalar dTL   = T_liq - T_sol;

        scalar H_old;
        if (T_curr <= T_sol)
            H_old = cp_p * T_curr;
        else if (T_curr >= T_liq)
            H_old = cp_p * T_curr + L_f;
        else
            H_old = cp_p * T_curr + L_f * (T_curr - T_sol) / dTL;

        // Eq.(10): explicit Euler step on enthalpy
        const scalar H_new = H_old + (Q_laser - Q_rad - Q_conv) * dt / m;

        if (H_new <= H_sol)
            T_curr = H_new / cp_p;
        else if (H_new >= H_liq)
            T_curr = (H_new - L_f) / cp_p;
        else
            T_curr = (H_new * dTL + L_f * T_sol) / (cp_p * dTL + L_f);  // mushy zone [Eq.(11)]

        T_persistent.set(key, T_curr);
        ownedKeys[owned]  = key;
        ownedTemps[owned] = T_curr;
        owned++;

        if (T_curr > T_max_local)
        {
            T_max_local   = T_curr;
            pos_Tmax_local = p.position();
        }
        const scalar dT_step = mag(T_curr - T_old);
        if (dT_step > dT_max_local)
        {
            dT_max_local    = dT_step;
            pos_dTmax_local = p.position();
        }

        if (Q_laser > 1.0e-6) nHeated++;
        
        
        // Eqs.(12-13): explicit Euler stability check
        // tau = m*cp / beta_eff;  beta_eff = A_surf*(h + 4*eps*sigma*T^3)
        // The 4T^3 term is the linearised radiation coefficient dQ_rad/dT.
        const scalar beta_eff = A_surf * (h_conv + 4.0*epsilon_p*sigma_SB*pow3(T_curr));
        const scalar tau_thermal = m * cp_p / max(beta_eff, SMALL);
        if (dt > 2.0 * tau_thermal)
        {
            WarningInFunction
                << "Particle " << key << ": dt=" << dt
                << " > 2*tau=" << 2.0*tau_thermal
                << " — explicit Euler may be unstable!" << nl;
        }
    }
    ownedKeys.resize(owned);
    ownedTemps.resize(owned);

    // -----------------------------------------------------------------------
    // Step 2: MPI sync -- broadcast post-computation values from current owner.
    // No if-not-found guard: the current owner always wins, handles migration.
    // -----------------------------------------------------------------------
    List<labelList> allKeys(Pstream::nProcs());
    List<scalarList> allTemps(Pstream::nProcs());

    allKeys[Pstream::myProcNo()]  = ownedKeys;
    allTemps[Pstream::myProcNo()] = ownedTemps;

    Pstream::gatherList(allKeys);
    Pstream::gatherList(allTemps);
    Pstream::broadcastList(allKeys);
    Pstream::broadcastList(allTemps);

    forAll(allKeys, procI)
    {
        forAll(allKeys[procI], i)
        {
            T_persistent.set(allKeys[procI][i], allTemps[procI][i]);
        }
    }

    // -----------------------------------------------------------------------
    // Per-step statistics output
    // -----------------------------------------------------------------------
    {
        const label  nTotal_g  = returnReduce(owned,   sumOp<label>());
        const label  nHeated_g = returnReduce(nHeated, sumOp<label>());

        // Gather per-proc dT_max and corresponding position, then find global argmax
        scalarList  all_dT(Pstream::nProcs(), scalar(0));
        List<vector> all_pos(Pstream::nProcs(), vector::zero);
        all_dT[Pstream::myProcNo()]  = dT_max_local;
        all_pos[Pstream::myProcNo()] = pos_dTmax_local;
        Pstream::gatherList(all_dT);
        Pstream::gatherList(all_pos);
        Pstream::broadcastList(all_dT);
        Pstream::broadcastList(all_pos);

        scalar dT_max_g    = 0.0;
        vector pos_dTmax_g = vector::zero;
        forAll(all_dT, procI)
        {
            if (all_dT[procI] > dT_max_g)
            {
                dT_max_g    = all_dT[procI];
                pos_dTmax_g = all_pos[procI];
            }
        }

        // Gather per-proc T_max and corresponding position, then find global argmax
        // Use running-maximum (same pattern as dT_max) to avoid floating-point
        // comparison pitfalls with T_max_g - SMALL at large T values.
        scalarList   all_Tmax(Pstream::nProcs(), scalar(-GREAT));
        List<vector> all_posTmax(Pstream::nProcs(), vector::zero);
        all_Tmax[Pstream::myProcNo()]    = owned > 0 ? T_max_local : scalar(-GREAT);
        all_posTmax[Pstream::myProcNo()] = pos_Tmax_local;
        Pstream::gatherList(all_Tmax);
        Pstream::gatherList(all_posTmax);
        Pstream::broadcastList(all_Tmax);
        Pstream::broadcastList(all_posTmax);

        scalar T_max_g     = scalar(-GREAT);
        vector pos_Tmax_g  = vector::zero;
        forAll(all_Tmax, procI)
        {
            if (all_Tmax[procI] > T_max_g)
            {
                T_max_g    = all_Tmax[procI];
                pos_Tmax_g = all_posTmax[procI];
            }
        }

        if (nTotal_g > 0)
        {
            Info<< "particleInfo [t=" << runTime.timeName() << "]:"
                << "  N=" << nTotal_g
                << "  heated=" << nHeated_g
                << "  T_max=" << T_max_g << " K"
                << "  pos(T_max)=" << pos_Tmax_g
                << "  dT_max=" << dT_max_g << " K"
                << "  pos(dT_max)=" << pos_dTmax_g
                << endl;
        }
    }

    // -----------------------------------------------------------------------
    // Step 3: Write "T" Lagrangian field to disk at write intervals only
    // -----------------------------------------------------------------------
    if (runTime.writeTime())
    {
        const word cloudName = kinematicCloudName;
        const label nLocal   = kinematicCloud.size();

        IOField<scalar> T_field
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                "lagrangian"/cloudName,
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            nLocal
        );

        label pI = 0;
        forAllConstIter(basicKinematicTypeCloud, kinematicCloud, iter)
        {
            const auto& p = iter();
            const label key = p.origProc() * 1000000 + p.origId();
            const label mId = (p.typeId() >= 0 && mat_Tinit.found(p.typeId()))
                             ? p.typeId() : 0;
            T_field[pI++] = T_persistent.lookup(key, mat_Tinit[mId]);
        }

        T_field.write();

        const label nTotal = returnReduce(nLocal, sumOp<label>());
        const label nHeatTotal = returnReduce(nHeated, sumOp<label>());
        Info<< "particleInfo: " << nTotal << " particles, "
            << nHeatTotal << " heated at t=" << runTime.timeName() << endl;
    }
}
