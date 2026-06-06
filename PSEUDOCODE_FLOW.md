# LaserAttDPTDMPFoam 求解器伪代码流程
## OpenFOAM v10 激光-粒子-气体耦合求解器

---

## 1. 整体架构伪代码

```
PROGRAM LaserAttDPTDMPFoam
│
├─ SUBROUTINE initialization()
│  ├─ Read mesh
│  ├─ Read mesh motion controls
│  ├─ Read fields (U, p, alpha)
│  ├─ Create particle clouds
│  ├─ Create fv schemes and solvers
│  ├─ Initialize laser fields
│  ├─ Create continuity error tracker
│  └─ Print header messages
│
├─ WHILE pimple.run(runTime) DO
│  │
│  ├─ TIME STEP INITIALIZATION (readDyMControls, CourantNo, setDeltaT)
│  │  ├─ Read dynamic mesh controls
│  │  ├─ Compute Courant number
│  │  ├─ Adjust time step if needed
│  │  └─ runTime++  [t → t + Δt]
│  │
│  ├─ PHASE 1: MESH UPDATE (动态网格)
│  │  ├─ mesh.update()           [拓扑和映射]
│  │  ├─ clouds.storeGlobalPositions()  [保存全局粒子位置]
│  │  ├─ mesh.move()             [网格运动]
│  │  └─ IF mesh.changing() THEN
│  │     ├─ correctPhic()        [修正通量]
│  │     └─ meshCourantNo()      [网格Courant数]
│  │     END IF
│  │
│  ├─ PHASE 2: PROPERTY UPDATE (物性更新)
│  │  ├─ continuousPhaseViscosity.correct()
│  │  └─ muc = rhoc * nu
│  │
│  ├─ PHASE 3: PARTICLE EVOLUTION (粒子演化 - Lagrangian)
│  │  ├─ CALL particleCloudEvolution()  [详见第3节]
│  │  │  └─ 内部：更新速度、位置、温度、碰撞等
│  │  │
│  │  ├─ alphac = max(1.0 - clouds.theta(), alphacMin)
│  │  │  └─ [气体体积分数 ← 1 - 粒子体积分数]
│  │  │
│  │  ├─ alphac.correctBoundaryConditions()
│  │  ├─ alphacf = interpolate(alphac)     [面值]
│  │  └─ alphaPhic = alphacf * phic        [修正通量]
│  │
│  ├─ PHASE 4: CLOUD FORCE COMPUTATION (粒子作用力)
│  │  ├─ cloudSU = clouds.SU(Uc)  [矩阵形式]
│  │  │  └─ [计算所有粒子受气体曳力、压力梯度力、重力]
│  │  │
│  │  └─ CASE cloudForceSplit OF
│  │     ├─ faceExplicitCellImplicit: (默认)
│  │     │  ├─ cloudSUu = -cloudSU.source() / mesh.V()  [显式部分→面]
│  │     │  ├─ cloudSUp = 0                             [无隐式修正]
│  │     │  └─ cloudSU.source() = 0                     [隐式部分保留→单元]
│  │     │
│  │     ├─ faceExplicitCellLagged:
│  │     │  ├─ cloudSUu = (cloudSU.diag()*Uc - cloudSU.source()) / mesh.V()
│  │     │  ├─ cloudSUp = 0
│  │     │  └─ cloudSU.source() = cloudSU.diag() * Uc    [滞后修正]
│  │     │
│  │     └─ faceImplicit:
│  │        ├─ cloudSUu = -cloudSU.source() / mesh.V()
│  │        ├─ cloudSUp = cloudSU.diag() / mesh.V()      [隐式到面]
│  │        └─ cloudSU.diag() = 0; cloudSU.source() = 0
│  │     END CASE
│  │
│  ├─ PHASE 5: GASEOUS FLOW SOLVING (气体求解 - Eulerian)
│  │  ├─ WHILE pimple.loop() DO  [PIMPLE外循环]
│  │  │  │
│  │  │  ├─ fvModels.correct()
│  │  │  │
│  │  │  ├─ CALL momentumEquationSolver()  [详见第4节]
│  │  │  │  ├─ Build fvVectorMatrix UcEqn(...)
│  │  │  │  ├─ UcEqn.relax()
│  │  │  │  ├─ IF pimple.momentumPredictor() THEN
│  │  │  │  │  └─ solve(UcEqn == rhs)
│  │  │  │  │  END IF
│  │  │  │  └─ Uc.correctBoundaryConditions()
│  │  │  │
│  │  │  └─ WHILE pimple.correct() DO  [PISO内循环]
│  │  │     ├─ CALL pressureCorrectionSolver()  [详见第5节]
│  │  │     │  ├─ Build fvScalarMatrix pEqn(...)
│  │  │     │  ├─ pEqn.solve()
│  │  │     │  ├─ phic = phiHbyASp - pEqn.flux() / alphacf
│  │  │     │  ├─ p.relax()
│  │  │     │  └─ Uc = HbyA + correction
│  │  │     │  END IF
│  │  │     └─ Uc.correctBoundaryConditions()
│  │  │     END WHILE
│  │  │
│  │  │  └─ IF pimple.turbCorr() THEN
│  │  │     └─ continuousPhaseTurbulence.correct()
│  │  │     END IF
│  │  └─ END WHILE
│  │
│  ├─ PHASE 6: OUTPUT TO DISK (输出)
│  │  └─ runTime.write()
│  │     └─ [write: 0.001/U, p, alpha, lagrangian/cloud/positions, d, U, ...]
│  │
│  ├─ PHASE 7: LASER COMPUTATIONS (激光计算 - 仅在输出步)
│  │  ├─ IF runTime.outputTime() THEN
│  │  │  │
│  │  │  ├─ CALL laserAttenuationComputation()  [详见第6节]
│  │  │  │  ├─ Read LaserProperties dictionary
│  │  │  │  ├─ IF enableAttenuation THEN
│  │  │  │  │  ├─ Parallel communication: 收集全局粒子数据
│  │  │  │  │  ├─ Identify cells in laser beam
│  │  │  │  │  ├─ Sort cells along laser direction
│  │  │  │  │  ├─ FOR each cell in beam DO
│  │  │  │  │  │  └─ laserAttenuationFactor[cell] = ∏(exp(-σ_ext/A_beam))
│  │  │  │  │  │  END FOR
│  │  │  │  │  ├─ laserAttenuationFactor.correctBoundaryConditions()
│  │  │  │  │  └─ laserAttenuationFactor.write()
│  │  │  │  │  ELSE
│  │  │  │  │  ├─ laserAttenuationFactor = 1.0  [无衰减]
│  │  │  │  │  END IF
│  │  │  │
│  │  │  ├─ CALL laserIntensityComputation()   [详见第7节]
│  │  │  │  ├─ Read laser parameters (HS_a, HS_Q, gaussianK, etc.)
│  │  │  │  ├─ FOR each cell DO
│  │  │  │  │  ├─ IF cell has particles THEN
│  │  │  │  │  │  ├─ Calculate distance r to laser axis
│  │  │  │  │  │  ├─ IF r < HS_a THEN
│  │  │  │  │  │  │  ├─ Q_base = (HS_Q / πa²) * exp(-gaussianK * r²/a²)
│  │  │  │  │  │  │  ├─ laserIntensity[cell] = Q_base * laserAttenuationFactor[cell]
│  │  │  │  │  │  │  └─ Calculate energy deposition
│  │  │  │  │  │  │  ELSE
│  │  │  │  │  │  │  └─ laserIntensity[cell] = 0
│  │  │  │  │  │  │  END IF
│  │  │  │  │  │  END IF
│  │  │  │  │  END FOR
│  │  │  │  └─ laserIntensity.write()
│  │  │  │
│  │  │  └─ CALL particleTemperatureEvolution()  [详见第8节]
│  │  │     ├─ Read particle properties (sigma, T_amb, alpha_laser, etc.)
│  │  │     ├─ Calculate dt_total = t_now - t_last_write
│  │  │     ├─ Calculate nSubCycles = ceil(dt_total / dt_simulation)
│  │  │     │
│  │  │     ├─ FOR each particle p DO
│  │  │     │  ├─ Read previous temperature T_old[p] (temperature tracking)
│  │  │     │  ├─ FOR sub = 1 to nSubCycles DO
│  │  │     │  │  │
│  │  │     │  │  ├─ cellI = mesh.findCell(position[p])
│  │  │     │  │  ├─ Calculate distance r_to_axis = |p.position - laser_axis|
│  │  │     │  │  │
│  │  │     │  │  ├─ I_laser = (HS_Q / πa²) * exp(-gaussianK * r²/a²) * laserAttenuationFactor[cellI]
│  │  │     │  │  ├─ Q_laser = alpha_laser * I_laser * A_proj
│  │  │     │  │  ├─ Q_rad = epsilon * sigma * A_surf * (T^4 - T_amb^4)
│  │  │     │  │  ├─ Q_conv = h * A_surf * (T - T_fluid)  [if enabled]
│  │  │     │  │  │
│  │  │     │  │  ├─ dT = (Q_laser - Q_rad - Q_conv) / (m * cp) * dt_sub
│  │  │     │  │  └─ T_new = T_old + dT
│  │  │     │  │
│  │  │     │  │  END FOR
│  │  │     │  │
│  │  │     │  ├─ T_particle[p] = T_new
│  │  │     │  └─ OUTPUT(p_id, position, T_old, T_new, ΔT, Q_laser, Q_rad)
│  │  │     │
│  │  │     └─ END FOR
│  │  │     └─ T_particle.write()  [save to lagrangian/cloud/T_particle]
│  │  │
│  │  │  END IF  [outputTime]
│  │
│  ├─ PHASE 8: STATISTICS (性能统计)
│  │  ├─ Print continuity errors
│  │  ├─ Print execution time
│  │  └─ Print clock time
│  │
│  └─ END WHILE  [主时间循环]
│
└─ PRINT "End"

END PROGRAM
```

---

## 2. 初始化阶段详解

```
FUNCTION initialization()
│
├─ SET up mesh
│  ├─ Read "mesh" directory
│  ├─ mesh object ready for FV discretization
│  └─ Prepare mesh schemes and tools
│
├─ SET up fields (从0/目录读取)
│  ├─ U.air        [气体速度场]
│  ├─ p            [压力场]
│  ├─ phi.air      [通量场]
│  └─ alpha.air    [气体体积分数]
│
├─ CREATE continuous phase properties
│  ├─ Read physicalProperties (重要！)
│  │  └─ continuousPhaseName = "air"
│  │
│  ├─ rhoc = ρ_continuous  [气体密度]
│  ├─ muc = μ_continuous   [气体粘度]
│  └─ Create viscosity model (Newtonian/Non-Newtonian)
│
├─ CREATE particle clouds
│  ├─ CALL parcelCloudList(rhoc, Uc, muc, g)
│  │  └─ Reads constant/cloudProperties
│  │     ├─ Cloud type (collidingCloud, ...)
│  │     ├─ Particle forces models (WenYuDrag, gravity, ...)
│  │     ├─ Collision model (pairCollision, ...)
│  │     └─ Create numerical integrator
│  │
│  ├─ Set Alpha_c_min (minimum air volume fraction)
│  │  └─ From fvSolution: alpha.air → max → alphacMin = 1 - max
│  │
│  └─ Initialize alphac from cloud
│     └─ alphac = max(1.0 - clouds.theta(), alphacMin)
│
├─ CREATE turbulence model
│  ├─ Read turbulenceProperties
│  ├─ Create momentumTransportModel
│  └─ Initialize nu_t (eddy viscosity)
│
├─ CREATE pressure reference
│  ├─ Set reference cell
│  └─ Set reference pressure value
│
├─ CREATE laser fields
│  ├─ laserAttenuationFactor[0..1]  [初始化为1.0]
│  ├─ laserIntensity[W/m²]          [初始化为0]
│  ├─ Deposition[W/m²]              [初始化为0]
│  └─ LaserPower[W]                 [初始化为0]
│
├─ SET up PIMPLE control
│  ├─ nOuterCorrectors = 1
│  ├─ nCorrectors = 3 (PISO sub-steps)
│  └─ cloudForceSplit = faceExplicitCellImplicit
│
├─ INITIALIZE continuity error tracking
│  ├─ cumulativeContErr = 0.0
│  └─ Create continuity error field
│
└─ PRINT header messages
   └─ "Starting time loop"

END FUNCTION
```

---

## 3. 粒子云演化 (clouds.evolve())

```
FUNCTION particleCloudEvolution()
│
├─ FOR each parcel (time stepping)
│  │
│  ├─ READ parcel properties
│  │  ├─ position (x, y, z)
│  │  ├─ velocity (U_x, U_y, U_z)
│  │  ├─ diameter d [m]
│  │  ├─ density ρ_p [kg/m³]
│  │  ├─ mass m [kg]
│  │  ├─ temperature T [K] (若启用)
│  │  └─ cell location cellI
│  │
│  ├─ LOCATE particle in mesh
│  │  └─ cellI = mesh.findCell(position)
│  │     └─ IF not found THEN particle exits, mark for deletion
│  │     END IF
│  │
│  ├─ COMPUTE forces ON parcel
│  │  │
│  │  ├─ Drag force (Wen-Yu model for dense suspension)
│  │  │  ├─ U_rel = U_gas[cellI] - U_p
│  │  │  ├─ α_p = 粒子体积分数[cellI]
│  │  │  ├─ Re_p = ρ_gas * |U_rel| * d / μ_gas
│  │  │  │
│  │  │  ├─ IF Re_p < 1000 THEN
│  │  │  │  └─ C_d = (24/(Re_p*α_c)) * (1 + 0.15*Re_p^0.687)
│  │  │  │  ELSE
│  │  │  │  └─ C_d = 0.44
│  │  │  │  END IF
│  │  │  │
│  │  │  └─ F_drag = (3/4) * C_d * α_p * ρ_gas * |U_rel| * U_rel
│  │  │           * (α_p^α_p / α_c^2.5)
│  │  │
│  │  ├─ Gravity force
│  │  │  └─ F_gravity = m * g
│  │  │
│  │  ├─ Pressure gradient force
│  │  │  └─ F_pressure = -V_p * ∇p[cellI]
│  │  │
│  │  └─ Collision forces (if particle-particle collision detected)
│  │     └─ F_collision = contact forces based on collision model
│  │
│  ├─ UPDATE velocity (Lagrangian ODE)
│  │  └─ dU/dt = (F_drag + F_gravity + F_pressure + F_collision) / m
│  │     └─ U_new = U_old + dU/dt * Δt
│  │
│  ├─ UPDATE position
│  │  └─ dx/dt = U_p
│  │     └─ position_new = position_old + U_new * Δt
│  │
│  ├─ HANDLE collisions
│  │  ├─ FOR each neighboring parcel DO
│  │  │  ├─ COMPUTE distance between parcels
│  │  │  ├─ IF distance < (r_i + r_j) THEN
│  │  │  │  ├─ COMPUTE collision normal
│  │  │  │  ├─ APPLY normal impulse (coefficient e_n)
│  │  │  │  ├─ APPLY tangential impulse (coefficient μ_t)
│  │  │  │  └─ Update velocities both parcels
│  │  │  │  END IF
│  │  │  END FOR
│  │  │
│  │  └─ Handle wall collisions
│  │     ├─ IF particle near boundary THEN
│  │     │  ├─ COMPUTE distance to wall
│  │     │  ├─ APPLY wall rebound (coefficient e_wall)
│  │     │  └─ Update velocity
│  │     │  END IF
│  │
│  ├─ PROJECT onto cells (cell finding)
│  │  └─ Update cellI for next step
│  │
│  └─ INCREMENT parcel integration counter
│     └─ nParcels_integrated++
│
├─ COMPUTE cloud volume fraction field
│  ├─ clouds.theta()
│  │  └─ [mesh grid: α_p[cell] = Σ(V_parcel in cell) / V_cell]
│  │
│  └─ [Used later for alphac update]
│
├─ COMPUTE cloud source terms (SU matrix)
│  └─ clouds.SU(Uc)
│     └─ [For each cell: source = Σ(F_drag_on_parcels_in_cell)]
│        └─ Matrix form: f_source (explicit) + f_diag * U (implicit)
│
└─ END FUNCTION

Output:
  ├─ Updated parcel positions
  ├─ Updated parcel velocities
│  ├─ clouds.theta() field         (粒子体积分数)
│  ├─ clouds.SU() matrix           (源项)
│  └─ Parcel count (nParcels)
```

---

## 4. 动量方程求解 (UcEqn.H)

```
FUNCTION momentumEquationSolver()
│
├─ BUILD momentum equation matrix
│  │
│  │  fvVectorMatrix UcEqn(
│  │      fvm::ddt(alphac, Uc)
│  │    + fvm::div(alphaPhic, Uc)
│  │    - fvm::Sp(fvc::ddt(alphac) + fvc::div(alphaPhic), Uc)
│  │    + turbulence.divDevTau(Uc)
│  │   ==
│  │      (1.0/rhoc) * cloudSU
│  │  )
│  │
│  ├─ First term: ∂(α_c * ρ_c * U) / ∂t  [瞬态项]
│  ├─ Second term: ∇·(α_c * ρ_c * U ⊗ U)  [对流项]
│  ├─ Third term: -∇·(α_c * ρ_c * U) * U   [连续性修正]
│  ├─ Fourth term: ∇·(α_c * τ_c)           [粘性项]
│  └─ RHS: (1/ρ_c) * F_particle            [粒子作用力源项]
│
├─ APPLY relaxation
│  └─ UcEqn.relax(rho=0.7, U=0.3)  [稳定性因子]
│
├─ APPLY constraints
│  └─ fvConstraints.constrain(UcEqn)  [边界条件约束]
│
├─ COMPUTE inverse of diagonal
│  ├─ rAUc = 1.0 / UcEqn.A()                 [单纯逆]
│  ├─ rASpUc = 1.0 / (UcEqn.A() - cloudSUp*rhoc)  [改进逆，含隐式修正]
│  └─ rASpUcf = fvc::interpolate(rASpUc)     [插值到面]
│
├─ COMPUTE face flux contributions
│  ├─ phicSUSu = fvc::flux(rASpUc * cloudSUu / rhoc)  + ...  [显式部分]
│  └─ phicSUSp = fvc::interpolate(rASpUc * cloudSUp / rhoc)  [隐式部分]
│
├─ IF momentum predictor enabled THEN
│  │
│  └─ SOLVE momentum equation
│     └─ solve(
│        │  UcEqn
│        │ ==
│        │  fvc::reconstruct((phicSUSu + phicSUSp*phic)/rASpUcf - ∇p*mesh.magSf())
│        │ + (1/rhoc) * (fvm::Sp(cloudSUp, Uc) - cloudSUp*Uc)
│        └─ )
│        └─ [RESULT: Uc预测值]
│
│  ├─ APPLY constraints
│  │  └─ fvConstraints.constrain(Uc)
│  │
│  └─ CORRECT boundary conditions
│     └─ Uc.correctBoundaryConditions()
│
│  END IF
│
└─ END FUNCTION

Output:
  ├─ Updated velocity field Uc (preliminary)
  └─ Intermediate matrix elements for pressure correction
```

---

## 5. 压力修正方程求解 (pEqn.H)

```
FUNCTION pressureCorrectionSolver()
│
├─ COMPUTE H by A field (momentum residual)
│  ├─ HbyA = constrainHbyA(rAUc * UcEqn.H(), Uc, p)
│  │  └─ [H() = RHS - diagonal*U from matrix equation]
│  │
│  └─ HbyASp = rASpUc / rAUc * HbyA  [改进版本]
│
├─ BUILD corrected mass flux
│  ├─ phiHbyASp = fvc::flux(HbyASp) + alphacf*rASpUcf*fvc::ddtCorr(Uc, phic, Ucf)
│  │
│  └─ [Includes time derivative correction for BDF schemes]
│
├─ ADJUST for mesh motion (if any)
│  ├─ IF pressureRef.needReference() THEN
│  │  ├─ fvc::makeRelative(phiHbyASp, Uc)
│  │  ├─ adjustPhi(phiHbyASp, Uc, p)
│  │  └─ fvc::makeAbsolute(phiHbyASp, Uc)
│  │  END IF
│
├─ ADD particle source contribution to flux
│  └─ phiHbyASp += phicSUSu
│
├─ ENFORCE pressure boundary conditions
│  └─ constrainPressure(p, Uc, phiHbyASp, rASpUcf)
│
├─ BUILD pressure equation
│  │
│  │  fvScalarMatrix pEqn(
│  │      fvm::laplacian(alphacf * rASpUcf, p)
│  │     ==
│  │      fvc::ddt(alphac) + fvc::div(alphacf * phiHbyASp)
│  │  )
│  │
│  ├─ LHS: ∇·(α_c/ρ_c * ∇p)  [Poisson operator]
│  └─ RHS: ∂α_c/∂t + ∇·(α_c * flux)  [continuity source]
│
├─ SET reference pressure
│  └─ pEqn.setReference(pressureRef.refCell(), pressureRef.refValue())
│
├─ SOLVE pressure equation
│  └─ pEqn.solve()
│     └─ [RESULT: updated pressure p]
│
├─ IF final non-orthogonal iteration THEN
│  │
│  ├─ CORRECT mass flux
│  │  └─ phic = phiHbyASp - pEqn.flux() / alphacf
│  │
│  ├─ RELAX pressure
│  │  └─ p.relax()  [rho=relTol_value]
│  │
│  ├─ CORRECT velocity
│  │  ├─ Uc = HbyA + rAUc * fvc::reconstruct((phicSUSu + phicSUSp*phic - pEqn.flux()/alphacf) / rASpUcf)
│  │  ├─ Uc.correctBoundaryConditions()
│  │  └─ fvConstraints.constrain(Uc)
│  │
│  ├─ CORRECT face velocity (if moving mesh)
│  │  └─ fvc::correctUf(Ucf, Uc, phic)
│  │
│  └─ MAKE flux relative to mesh motion
│     └─ fvc::makeRelative(phic, Uc)
│  │
│  END IF
│
└─ END FUNCTION

Output:
  ├─ Updated pressure field p
  ├─ Corrected velocity field Uc
  └─ Corrected mass flux phic
```

---

## 6. 激光衰减计算 (LaserAttenuation.H)

```
FUNCTION laserAttenuationComputation()
│
├─ IF NOT runTime.outputTime() THEN
│  └─ RETURN (激光计算仅在输出步执行)
│  END IF
│
├─ READ laser properties
│  ├─ enableAttenuation: BOOL
│  ├─ extinctionEfficiency: scalar (Qext)
│  ├─ V_incident: vector
│  ├─ laserInitialPosition: vector
│  ├─ laserVelocity: vector
│  └─ HS_a: scalar (beam radius)
│
├─ IF NOT enableAttenuation THEN
│  ├─ laserAttenuationFactor = 1.0  (everywhere)
│  └─ RETURN
│  END IF
│
├─ READ particle cloud
│  ├─ CALL passiveParticleCloud(mesh, "cloud")
│  ├─ nParticles_local = particles.size()
│  └─ readField(particles, "d")  [diameter]
│
├─ PARALLEL COMMUNICATION (collect particles from all processors)
│  │
│  ├─ Extract local particles
│  │  ├─ localParticlePos[i] = particles[i].position()
│  │  └─ localParticleDiam[i] = diameter[i]
│  │
│  ├─ Gather from all processes
│  │  ├─ allProcParticlePos[myProcNo] = localParticlePos
│  │  ├─ Pstream.gatherList(allProcParticlePos)      [data to master]
│  │  ├─ Pstream.scatterList(allProcParticlePos)     [distribute back]
│  │  └─ Now all processes have complete particle list
│  │
│  └─ Flatten into global arrays
│     ├─ FOR each processor DO
│     │  ├─ globalParticlePos += allProcParticlePos[proc]
│     │  └─ globalParticleDiam += allProcParticleDiam[proc]
│     │  END FOR
│     └─ nParticles_total = globalParticlePos.size()
│
├─ BUILD cell-to-particles mapping (global)
│  └─ FOR each global particle p DO
│     ├─ cellI = mesh.findCell(globalParticlePos[p])
│     ├─ IF cellI valid THEN
│     │  └─ cellParticles[cellI].append(p)
│     │  END IF
│     END FOR
│
├─ CALCULATE laser position at current time
│  ├─ laserPos = laserInitPos + laserVel * runTime.value()
│  ├─ laserDirection = V_incident / |V_incident|
│  └─ beamRadius = HS_a
│
├─ IDENTIFY cells in laser beam
│  │
│  └─ FOR each cell DO
│     ├─ cellCenter = mesh.C()[cell]
│     ├─ relPos = cellCenter - laserPos
│     ├─ distToAxis = |relPos - (relPos·laserDir)*laserDir|   [perpendicular distance]
│     │
│     ├─ IF distToAxis < 2.0 * beamRadius THEN
│     │  ├─ beamCells.append(cell)
│     │  ├─ beamCellProj.append(cellCenter·laserDir)  [projection along beam]
│     │  END IF
│     END FOR
│
├─ SORT cells along laser direction (upstream → downstream)
│  └─ sortIndices = sortedOrder(beamCellProj)
│     └─ [Now cells ordered by position along laser beam]
│
├─ CALCULATE total extinction cross-section for each cell
│  │
│  └─ FOR each cell DO
│     ├─ crossSection[cell] = 0.0
│     ├─ FOR each particle p in cell DO
│     │  ├─ d = globalParticleDiam[p]
│     │  ├─ σ_single = Qext * π * d² / 4
│     │  └─ crossSection[cell] += σ_single
│     │  END FOR
│     END FOR
│
├─ COMPUTE cumulative attenuation (Layered approach)
│  │
│  ├─ beamArea = π * beamRadius²
│  │
│  └─ FOR i = 0 to beamCells.size()-1 DO
│     │  
│     ├─ cellI = beamCells[sortIndices[i]]
│     ├─ cumulativeAttenuation = 1.0  [start with no attenuation]
│     │
│     ├─ FOR j = 0 to i-1 DO  [loop over upstream cells]
│     │  │
│     │  ├─ upCellI = beamCells[sortIndices[j]]
│     │  │
│     │  ├─ Check if cells are in same radial "tube"
│     │  │  ├─ distToAxis_i = distance(cellI to laser axis)
│     │  │  ├─ distToAxis_j = distance(upCellI to laser axis)
│     │  │  ├─ radialDiff = |distToAxis_i - distToAxis_j|
│     │  │  └─ IF radialDiff < cellSize THEN   [in same tube]
│     │  │     │
│     │  │     ├─ Apply Beer-Lambert attenuation
│     │  │     ├─ upAttenuation = exp(-crossSection[upCellI] / beamArea)
│     │  │     └─ cumulativeAttenuation *= upAttenuation
│     │  │     │
│     │  │     END IF
│     │  │
│     │  END FOR
│     │
│     ├─ laserAttenuationFactor[cellI] = cumulativeAttenuation
│     │  └─ [Range: [0, 1], where 1=no attenuation, 0=complete blocking]
│     │
│     END FOR
│
├─ APPLY boundary conditions
│  └─ laserAttenuationFactor.correctBoundaryConditions()
│
└─ END FUNCTION

Output:
  └─ Updated field: laserAttenuationFactor[cell] ∈ [0, 1]
```

---

## 7. 激光强度计算 (LaserHS.H)

```
FUNCTION laserIntensityComputation()
│
├─ RESET intensity fields
│  ├─ Deposition *= 0.0     [已沉积能量复位]
│  ├─ LaserPower *= 0.0     [总功率复位]
│  └─ laserIntensity *= 0.0 [强度复位]
│
├─ CALCULATE cell geometric dimensions
│  │
│  └─ FOR each cell DO
│     ├─ Get cell vertices (points)
│     ├─ xDim[cell] = max(x) - min(x)  [X方向长度]
│     ├─ [可选: yDim, zDim]
│     END FOR
│
├─ READ laser parameters
│  ├─ HS_a: beam radius [m]
│  ├─ HS_Q: laser power [W]
│  ├─ gaussianK: gaussian parameter [1~3]
│  ├─ laserInitPos, laserVel, V_incident: geometry
│  ├─ HS_lg, HS_bg, HS_velocity: Y-direction scanning
│  └─ pi = π
│
├─ NORMALIZE laser direction
│  └─ laserDirection = V_incident / |V_incident|
│
├─ FOR each cell DO
│  │
│  ├─ IF cell has particles (alphac < 1.0) THEN
│  │  │
│  │  ├─ TRANSFORM cell center to local coordinate system
│  │  │  ├─ x_local = cell.x
│  │  │  ├─ y_local = cell.y - HS_lg + (HS_velocity * time)  [Y扫描]
│  │  │  └─ z_local = cell.z - HS_bg
│  │  │
│  │  ├─ CALCULATE distance to laser axis
│  │  │  └─ r = |(V1) × laserDirection|  [perpendicular distance]
│  │  │     where V1 = [x_local, y_local, z_local]
│  │  │
│  │  ├─ IF r < HS_a THEN  [in laser beam]
│  │  │  │
│  │  │  ├─ CALCULATE gaussian intensity distribution
│  │  │  │  ├─ Q_base = (2 * HS_Q) / (π * HS_a²) 
│  │  │  │  ├─ Q_base *= exp(-gaussianK * r² / HS_a²)
│  │  │  │  └─ [Result: [W/m²] at position r]
│  │  │  │
│  │  │  ├─ APPLY laser attenuation factor
│  │  │  │  ├─ Q_attenuated = Q_base * laserAttenuationFactor[cell]
│  │  │  │  └─ [Accounts for particle blocking upstream]
│  │  │  │
│  │  │  ├─ STORE laser intensity
│  │  │  │  └─ laserIntensity[cell] = Q_attenuated
│  │  │  │
│  │  │  ├─ CALCULATE energy deposition
│  │  │  │  ├─ Deposition[cell] += Q_attenuated / xDim[cell]
│  │  │  │  └─ LaserPower[cell] += Q_attenuated * mesh.V()[cell] / xDim[cell]
│  │  │  │
│  │  │  ELSE  [r >= HS_a, outside beam]
│  │  │  │
│  │  │  └─ laserIntensity[cell] = 0
│  │  │
│  │  │  END IF
│  │  │
│  │  ELSE  [no particles]
│  │  │
│  │  └─ laserIntensity[cell] = 0
│  │
│  │  END IF
│
│  ELSE  [cell.alphac >= 1.0, only gas]
│  │
│  └─ laserIntensity[cell] = 0
│
│  END IF
│
│  END FOR
│
├─ ACCUMULATE global statistics
│  ├─ TotalQ = gSum(LaserPower)     [total power deposited]
│  ├─ TotalIntensity = gSum(laserIntensity)  [total intensity]
│  └─ Print statistics
│
└─ END FUNCTION

Output:
  ├─ Updated field: laserIntensity[cell] [W/m²]
  ├─ Updated field: Deposition[cell]
  └─ Updated field: LaserPower[cell]
```

---

## 8. 粒子温度演化 (particleInfo.H)

```
FUNCTION particleTemperatureEvolution()
│
├─ IF NOT runTime.writeTime() THEN
│  └─ RETURN (仅在输出步计算粒子温度)
│  END IF
│
├─ READ static particle properties (cached after first call)
│  ├─ sigma = 5.67e-8 [Stefan-Boltzmann, W/(m²·K⁴)]
│  ├─ T_amb = 300.0   [ambient temperature, K]
│  ├─ alpha_laser = 0.4  [laser absorptivity]
│  ├─ epsilon = 0.3    [emissivity]
│  ├─ cp = 385.0       [specific heat, J/(kg·K)]
│  └─ rhoP = 8960.0    [particle density, kg/m³]
│
├─ READ dynamic laser parameters
│  ├─ HS_Q, HS_a, gaussianK
│  ├─ laserInitPos, laserVel, V_incident
│  └─ laserDirection = V_incident / |V_incident|
│
├─ CALCULATE time interval since last write
│  ├─ dt_total = runTime.value() - lastWriteTime
│  ├─ lastWriteTime = runTime.value()
│  ├─ IF dt_total < ε THEN dt_total = runTime.deltaTValue() [fallback]
│  │  END IF
│  │
│  ├─ Calculate number of sub-cycles
│  │  ├─ dt_sub = runTime.deltaTValue()
│  │  ├─ nSubCycles = max(1, ceil(dt_total / dt_sub))
│  │  └─ dt_sub = dt_total / nSubCycles  [actual sub time step]
│  │
│  └─ Info: Print dt information
│
├─ READ particle cloud
│  └─ CALL passiveParticleCloud(mesh, "cloud")
│
├─ FOR each particle p DO  [iterate cloud]
│  │
│  ├─ EXTRACT particle properties
│  │  ├─ cellI = p.cell()
│  │  ├─ globalPos = p.position()  [AUTO CONVERSION to global xyz]
│  │  ├─ diameter[p]
│  │  └─ origId[p], origProcId[p]  [for temperature tracking]
│  │
│  ├─ CALCULATE geometric properties
│  │  ├─ A_proj = π * (d/2)² = π*d²/4        [projected area]
│  │  ├─ A_surf = π * d²                     [surface area]
│  │  ├─ V = π/6 * d³                        [volume]
│  │  └─ m = rhoP * V                        [mass]
│  │
│  ├─ READ previous temperature (temperature tracking)
│  │  ├─ IF previous T_particle file exists THEN
│  │  │  ├─ Create particleKey = origProcId * 1e6 + origId
│  │  │  ├─ T_map[particleKey] = T_prev_value
│  │  │  ├─ IF tracked particle THEN
│  │  │  │  └─ T_curr = T_map[particleKey]   [inherit history]
│  │  │  │  ELSE
│  │  │  │  └─ T_curr = T_amb                [new particle]
│  │  │  │  END IF
│  │  │  │
│  │  │  ELSE  [no previous file]
│  │  │  │
│  │  │  └─ T_curr = T_amb  [first time]
│  │  │
│  │  │  END IF
│  │
│  ├─ STORE initial temperature
│  │  └─ T_initial = T_curr
│  │
│  ├─ GET attenuation factor
│  │  ├─ cellI = mesh.findCell(globalPos)
│  │  ├─ IF cellI >= 0 THEN
│  │  │  └─ attenuationFactor = laserAttenuationFactor[cellI]
│  │  │  ELSE
│  │  │  └─ attenuationFactor = 1.0  [outside domain]
│  │  │  END IF
│  │
│  ├─ SUB-CYCLE integration (improved accuracy)  ▼▼▼
│  │  │
│  │  └─ FOR sub = 1 to nSubCycles DO
│  │     │
│  │     ├─ CALCULATE distance to laser axis
│  │     │  ├─ laserPos = laserInitPos + laserVel * runTime.value()
│  │     │  ├─ relPos = globalPos - laserPos
│  │     │  └─ r = |relPos × laserDirection|  [perpendicular distance]
│  │     │
│  │     ├─ CALCULATE laser intensity at particle location
│  │     │  ├─ I_base = (HS_Q / (π * HS_a²)) * exp(-gaussianK * r² / HS_a²)
│  │     │  ├─ I_laser = I_base * attenuationFactor
│  │     │  └─ [Result: W/m²]
│  │     │
│  │     ├─ CALCULATE heat inputs/outputs
│  │     │  │
│  │     │  ├─ Energy input from laser:
│  │     │  │  └─ Q_laser = alpha_laser * I_laser * A_proj  [W]
│  │     │  │
│  │     │  ├─ Radiative cooling (Stefan-Boltzmann):
│  │     │  │  └─ Q_rad = epsilon * sigma * A_surf * (T_curr⁴ - T_amb⁴)  [W]
│  │     │  │
│  │     │  └─ Convective cooling (Newton's law):
│  │     │     └─ Q_conv = h * A_surf * (T_curr - T_fluid[cellI])  [optional]
│  │     │         where h = convection coefficient [W/(m²·K)]
│  │     │
│  │     ├─ CALCULATE net heat and temperature change
│  │     │  ├─ Q_net = Q_laser - Q_rad - Q_conv  [W]
│  │     │  ├─ dT = (Q_net * dt_sub) / (m * cp)  [K]
│  │     │  └─ T_curr = T_curr + dT
│  │     │
│  │     └─ END FOR [sub-cycle]
│  │
│  ├─ STORE final temperature
│  │  └─ T_particle[p] = T_curr
│  │
│  ├─ CALCULATE temperature rise
│  │  └─ temp_rise = T_curr - T_initial
│  │
│  └─ OUTPUT detailed heating information (if heating significant)
│     └─ Info: "Particle_ID, position, T_initial, T_final, ΔT, Q_laser, Q_rad"
│
│  END FOR  [particle loop]
│
├─ WRITE temperature field to disk
│  └─ T_particle.write()
│     └─ [Saved to: <time>/lagrangian/cloud/T_particle]
│
├─ UPDATE tracking variables
│  ├─ lastWrittenTime = runTime.timeName()
│  └─ hasLastWrittenTime = true
│
└─ END FUNCTION

Output:
  ├─ Updated field: T_particle[parcel]  [K]
  │  └─ Saved to lagrangian/cloud/T_particle
  └─ Details: heating statistics for each parcel
```

---

## 9. 完整时间步流程总结

```
Main Time Loop Iteration:
═══════════════════════════════════════════════════════════════════

STEP 1: Time Control (< 1% of time)
  ├─ readDyMControls()
  ├─ CourantNo()          [Check CFL]
  ├─ setDeltaT()          [Adjust Δt if needed]
  └─ runTime++            [t = t + Δt]

STEP 2: Mesh Update (< 2% of time)
  ├─ mesh.update()
  ├─ clouds.storeGlobalPositions()
  └─ mesh.move()

STEP 3: Property Update (< 1% of time)
  ├─ continuousPhaseViscosity.correct()
  └─ muc = rhoc * nu

STEP 4: Particle Evolution (Lagrangian) (≈ 15-25% of time)
  └─ clouds.evolve()      [Update U_p, x_p, handle collisions]

STEP 5: Volume Fraction Update (< 1% of time)
  ├─ alphac = max(1 - clouds.theta(), alphacMin)
  ├─ alphacf = interpolate(alphac)
  └─ alphaPhic = alphacf * phic

STEP 6: Cloud Force Computation (< 1% of time)
  ├─ cloudSU = clouds.SU(Uc)
  └─ Switch/distribute sources (faceExplicitCellImplicit, etc.)

STEP 7: Gas Flow Solving (≈ 50-70% of time) [PIMPLE loop]
  │
  ├─ WHILE pimple.loop() DO  [typically 1-3 outer iterations]
  │  │
  │  ├─ Momentum Equation (UcEqn.H)
  │  │  ├─ Build matrix with Uc, alphac, cloudSU
  │  │  ├─ IF momentumPredictor THEN solve
  │  │  └─ Calculate residuals
  │  │
  │  └─ WHILE pimple.correct() DO  [typically 2-3 inner iterations]
  │     │
  │     ├─ Pressure Correction (pEqn.H)
  │     │  ├─ Build Poisson equation for p'
  │     │  ├─ Solve for pressure corrections
  │     │  └─ Update Uc and phic
  │     │
  │     └─ Turbulence Correction
  │        └─ turbulence.correct()
  │
  └─ END WHILE (achieved convergence criteria)

STEP 8: Output to Disk (< 2% of time, only @ writeInterval)
  └─ runTime.write()
     └─ Fields: U.air, p, alpha.air, phi.air, laserIntensity, ...
     └─ Lagrangian: positions, d, U, T_particle, origId, ...

STEP 9: Laser Computations (≈ 5-10% of time, only @ writeInterval)
  │
  ├─ LaserAttenuation.H
  │  ├─ Parallel communication (gather particles)
  │  ├─ Build cell-to-particle mapping
  │  ├─ Calculate extinction cross-section
  │  └─ Compute layered attenuation: f_att[cell]
  │
  ├─ LaserHS.H
  │  ├─ FOR each cell
  │  │  ├─ Calculate distance r to laser axis
  │  │  ├─ Q_intensity = gaussian(r) * f_att
  │  │  └─ Store laserIntensity
  │  END FOR
  │
  └─ particleInfo.H
     ├─ FOR each particle
     │  ├─ FOR sub = 1 to nSubCycles
     │  │  ├─ Q_laser = alpha_laser * I * A_proj
     │  │  ├─ Q_rad = epsilon * sigma * (T⁴ - T_amb⁴)
     │  │  ├─ dT = ΔT / (m * cp)
     │  │  └─ T_new = T_old + dT
     │  END FOR
     │  └─ T_particle[p] = T_new
     END FOR
     └─ Save T_particle to disk

STEP 10: Statistics (< 1% of time)
  ├─ Print continuity errors
  ├─ Print execution times
  └─ END OF TIME STEP

Total Time per step: ~100% (proportions approximate, problem-dependent)
═══════════════════════════════════════════════════════════════════
```

---

## 10. 关键数据流图

```
时间层 t^n₊₁ 的结构：

入口:
  ├─ U^n, p^n, α_c^n         [来自t^n]
  ├─ U_p^n, x_p^n, T_p^n     [粒子数据 t^n]
  └─ Δt                      [时间步长]

并行处理（核心耦合）:
  │
  ├─ BRANCH 1 (Lagrangian)
  │  └─ particles.evolve()
  │     ├─ α_p^{n+1} (clouds.theta)
  │     ├─ U_p^{n+1}
  │     ├─ x_p^{n+1}
  │     └─ cloudSU    [⭐ 关键源项]
  │
  ├─ BRANCH 2 (Gas phase α update)
  │  ├─ α_c^{n+1} = 1 - α_p^{n+1}
  │  ├─ α_c^f (interpolation to faces)
  │  └─ α_c·φ (flux modification)
  │
  └─ BRANCH 3 (Gas phase solve)
     ├─ U^{n+1}        [from UcEqn + pEqn]
     ├─ p^{n+1}
     └─ φ^{n+1}

顺序输出（仅writeInterval）:
  │
  ├─ runTime.write()  [disk I/O]
  │  └─ Read: positions from 0.001/lagrangian/cloud/
  │
  ├─ LaserAttenuation.H  [Post-processing]
  │  ├─ Read: global particles
  │  ├─ Compute: f_att[cell]
  │  └─ Write: laserAttenuationFactor
  │
  ├─ LaserHS.H  [Post-processing]
  │  ├─ Compute: Q_laser[cell]
  │  └─ Write: laserIntensity
  │
  └─ particleInfo.H  [Post-processing]
     ├─ Compute: T_p^{n+1}
     └─ Write: T_particle

输出:
  ├─ All fields at t^{n+1}
  ├─ Ready for next iteration or restart
  └─ Post-processing results available
```

---

## 11. 优化与并行策略

```
PARALLELIZATION STRATEGY:

Domain Decomposition:
  ├─ Mesh partitioned into domains (one per processor)
  ├─ Particles distributed to processes with local mesh regions
  └─ Ghost layer communication for boundary cells

时间步中的并行部分:
  │
  ├─ ✓ clouds.evolve() [local only]
  │
  ├─ ✓ UcEqn/pEqn solve [local + MPI communication]
  │
  ├─ ✓ turbulence operations [local only]
  │
  └─ ✗ LaserAttenuation [requires GLOBAL particle data - Pstream.gatherList]

通信开销:
  ├─ Minimal during time stepping
  ├─ One Pstream.gatherList at each writeInterval (expensive)
  ├─ Solution: Enable attenuation only @ writeInterval
  └─ Alternative: Approximate local attenuation (not implemented)

性能瓶颈:
  ├─ Gas solve (≈60% of CPU time) - unavoidable
  ├─ Particle evolution (≈20% of CPU time) - decent scaling
  ├─ Laser computation (≈10% at writeInterval) - communication bound
  └─ I/O (≈10% at writeInterval) - file system limited

RECOMMENDED SETTINGS:
  ├─ nProcs: 8-64 (depending on mesh size)
  ├─ Particles per domain: 1k-10k (optimal)
  ├─ writeInterval: 0.005-0.01 s (to minimize laser overhead)
  └─ nOuterCorrectors: 1 (usually sufficient)
```

---

## 12. 故障排除指南

```
常见问题诊断:

问题 1: 数值不稳定 (NaN, Inf)
├─ 检查项:
│  ├─ Courant数: max (|U|·Δt/h) < 1.0  ✓
│  ├─ 时间步: try smaller Δt
│  ├─ 粒子体积分数: α_p < 0.6 ✓
│  ├─ 曳力: ✓ check Wen-Yu model at high α
│  └─ 激光功率: ✗ too high → excessive heating → divergence
└─ 解决:
   ├─ Reduce maxCo
   ├─ Check laser parameters (HS_Q should be physical)
   └─ Enable writeInterval output for debugging

问题 2: 粒子温度异常
├─ 检查项:
│  ├─ enableAttenuation = true? ✓
│  ├─ 衰减因子范围? [0, 1] ✓
│  ├─ 激光强度 > 0? Check laserIntensity field
│  └─ dt_sub > SMALL? Sub-cycling active?
└─ 解决:
   ├─ Verify LaserAttenuation output
   ├─ Check particle.pos vs laser.pos distance
   └─ Increase nSubCycles if temp oscillates

问题 3: 模拟速度慢
├─ 可能原因:
│  ├─ Mesh too fine (too many cells/particles)
│  ├─ PIMPLE nOuterCorrectors too high
│  ├─ Pressure tolerance too tight (1e-7 → 1e-6)
│  └─ Laser computation at each timestep (should be @ writeInterval only)
└─ 优化:
   ├─ Reduce mesh density if acceptable
   ├─ Set nOuterCorrectors = 1
   ├─ Increase relTol values (0.01-0.1)
   └─ Verify LaserAttenuation only runs @ writeInterval

问题 4: MPI通信错误
├─ 症状: Parallel communication fails
└─ 原因/解决:
   ├─ Array sizes mismatch: check clouds.size() consistency
   ├─ Gatherlist/scatterlist order: verify mesh partitioning
   └─ Solution: Run serial test first (1 process)
```

---

## 总结

各模块在时间步中的执行顺序和相对计算量：

```
时间步内 CPU时间 分配 (大约):

┌────────────────────────────────────────────┐
│ Time Step Total ≈ 100%                    │
├────────────────────────────────────────────┤
│ Step 1-3: Initialization        ≈ 2%     │
│ Step 4: Particle Evolution      ≈ 20%    │  ← Lagrangian
│ Step 5-6: Volume/Force Update   ≈ 2%     │
│ Step 7: Gas Flow Solve          ≈ 60%    │  ← Eulerian (PIMPLE+Poisson)
│ Step 8: Disk I/O               ≈ 5%     │  ← Only @ writeInterval
│ Step 9: Laser Computation      ≈ 10%    │  ← Only @ writeInterval
│ Step 10: Statistics             ≈ 1%     │
└────────────────────────────────────────────┘

Note:
  • Percentages vary with particle count and mesh resolution
  • LaserComputation (Step 9) only every ~100-1000 timesteps
  • Gas solve (Step 7) is computation-bound, parallelizes well
  • Particle evolution (Step 4) scales linearly with particle count
```

---

**文档版本：1.0**  
**日期：2024-03-12**  
**适用于：LaserAttDPTDMPFoam (OpenFOAM v10)**
