/* This file is part of the Palabos library.
 *
 * The Palabos software is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* \file
 * External flow in a domain defined by an STL file, with a void acting as the obstacle.
 * Modified to use a single STL for the domain, where the void is the solid obstacle.
 * Features:
 * - Loads domain geometry from STL file.
 * - Voxelizes domain with void as solid (obstacle).
 * - Applies boundary conditions on STL outer surfaces and void.
 * - Supports off-lattice boundary conditions, sponge zones, and particle streamlines.
 * - Computes force on the void (obstacle).
 * */

 #include "palabos3D.h"
 #include "palabos3D.hh"
 
 using namespace plb;
 using namespace std;
 
 typedef double T;
 typedef Array<T, 3> Velocity;
 #define DESCRIPTOR descriptors::D3Q19Descriptor
 typedef DenseParticleField3D<T, DESCRIPTOR> ParticleFieldT;
 
 #define PADDING 8
 
 static std::string outputDir("./tmp/");
 
 // Structure to hold user-defined and derived parameters
 struct Param {
     T nu;                         // Kinematic viscosity
     T lx, ly, lz;                 // Size of computational domain (derived from STL)
     bool freeSlipWall;            // Free-slip on void surface (obstacle)?
     bool lateralFreeSlip;         // Free-slip lateral boundaries?
     T maxT, statT, imageT, vtkT;  // Event times (physical units)
     plint resolution;             // Lattice nodes along reference length
     T inletVelocity;              // Inlet x-velocity (physical units)
     T uLB;                        // Velocity in lattice units
     bool useSmago;                // Use Smagorinsky LES model?
     T cSmago;                     // Smagorinsky parameter
     plint nx, ny, nz;             // Grid resolution
     T omega;                      // Relaxation parameter
     T dx, dt;                     // Discrete space and time steps
     plint maxIter, statIter;      // Event iterations
     plint imageIter, vtkIter;
     bool useParticles;            // Simulate particles?
     int particleTimeFactor;       // Particle integration time step factor
     T particleProbabilityPerCell; // Particle injection probability
     T cutOffSpeedSqr;             // Particle removal criterion
     int maxNumParticlesToWrite;   // Max particles in VTK output
     T outletSpongeZoneWidth;      // Outlet sponge zone width
     plint numOutletSpongeCells;   // Sponge zone lattice nodes
     int outletSpongeZoneType;     // Sponge zone type (0=Viscosity, 1=Smagorinsky)
     T targetSpongeCSmago;         // Target Smagorinsky parameter for sponge
     plint initialIter;            // Iterations for inlet velocity ramp-up
 
     // Boundary regions (adjusted for STL faces)
     Box3D inlet, outlet, lateral1, lateral2, lateral3, lateral4;
 
     std::string geometry_fname;
 
     Param() { }
 
     Param(std::string xmlFname)
     {
         XMLreader document(xmlFname);
         document["geometry"]["filename"].read(geometry_fname);
         document["geometry"]["freeSlipWall"].read(freeSlipWall);
         document["geometry"]["lateralFreeSlip"].read(lateralFreeSlip);
 
         // Derive domain size from STL bounding box
         TriangleSet<T> triangleSet(geometry_fname, DBL);
         Cuboid<T> bCuboid = triangleSet.getBoundingCuboid();
         lx = bCuboid.upperRightCorner[0] - bCuboid.lowerLeftCorner[0];
         ly = bCuboid.upperRightCorner[1] - bCuboid.lowerLeftCorner[1];
         lz = bCuboid.upperRightCorner[2] - bCuboid.lowerLeftCorner[2];
 
         document["numerics"]["nu"].read(nu);
         document["numerics"]["inletVelocity"].read(inletVelocity);
         document["numerics"]["resolution"].read(resolution);
         document["numerics"]["uLB"].read(uLB);
         document["numerics"]["useSmago"].read(useSmago);
         if (useSmago) {
             document["numerics"]["cSmago"].read(cSmago);
         }
 
         document["numerics"]["useParticles"].read(useParticles);
         if (useParticles) {
             document["numerics"]["particleTimeFactor"].read(particleTimeFactor);
             document["numerics"]["particleProbabilityPerCell"].read(particleProbabilityPerCell);
             document["numerics"]["cutOffSpeedSqr"].read(cutOffSpeedSqr);
             document["numerics"]["maxNumParticlesToWrite"].read(maxNumParticlesToWrite);
         }
 
         document["numerics"]["outletSpongeZoneWidth"].read(outletSpongeZoneWidth);
         std::string zoneType;
         document["numerics"]["outletSpongeZoneType"].read(zoneType);
         if ((util::tolower(zoneType)).compare("viscosity") == 0) {
             outletSpongeZoneType = 0;
         } else if ((util::tolower(zoneType)).compare("smagorinsky") == 0) {
             outletSpongeZoneType = 1;
         } else {
             pcout << "The sponge zone type must be either \"Viscosity\" or \"Smagorinsky\"." << std::endl;
             exit(-1);
         }
         document["numerics"]["targetSpongeCSmago"].read(targetSpongeCSmago);
 
         document["numerics"]["initialIter"].read(initialIter);
 
         document["output"]["maxT"].read(maxT);
         document["output"]["statT"].read(statT);
         document["output"]["imageT"].read(imageT);
         document["output"]["vtkT"].read(vtkT);
 
         computeLBparameters();
     }
 
     void computeLBparameters()
     {
         dx = ly / (resolution - 1.0); // Use ly as reference length
         dt = (uLB / inletVelocity) * dx;
         T nuLB = nu * dt / (dx * dx);
         omega = 1.0 / (DESCRIPTOR<T>::invCs2 * nuLB + 0.5);
         nx = util::roundToInt(lx / dx) + 1;
         ny = util::roundToInt(ly / dx) + 1;
         nz = util::roundToInt(lz / dx) + 1;
         maxIter = util::roundToInt(maxT / dt);
         statIter = util::roundToInt(statT / dt);
         imageIter = util::roundToInt(imageT / dt);
         vtkIter = util::roundToInt(vtkT / dt);
         numOutletSpongeCells = util::roundToInt(outletSpongeZoneWidth / dx);
 
         // Define boundary regions based on STL bounding box
         inlet = Box3D(0, 0, 0, ny - 1, 0, nz - 1); // x=0 face
         outlet = Box3D(nx - 1, nx - 1, 0, ny - 1, 0, nz - 1); // x=nx-1 face
         lateral1 = Box3D(1, nx - 2, 0, 0, 0, nz - 1); // y=0 face
         lateral2 = Box3D(1, nx - 2, ny - 1, ny - 1, 0, nz - 1); // y=ny-1 face
         lateral3 = Box3D(1, nx - 2, 1, ny - 2, 0, 0); // z=0 face
         lateral4 = Box3D(1, nx - 2, 1, ny - 2, nz - 1, nz - 1); // z=nz-1 face
     }
 
     Box3D boundingBox() const
     {
         return Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1);
     }
 
     T getInletVelocity(plint iIter)
     {
         static T pi = std::acos((T)-1.0);
         if (iIter >= initialIter) {
             return uLB;
         }
         if (iIter < 0) {
             iIter = 0;
         }
         return uLB * std::sin(pi * iIter / (2.0 * initialIter));
     }
 };
 
 Param param;
 
 // Apply boundary conditions to STL outer surfaces
 void outerDomainBoundaries(
     MultiBlockLattice3D<T, DESCRIPTOR> *lattice, MultiScalarField3D<T> *rhoBar,
     MultiTensorField3D<T, 3> *j, OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc)
 {
     Array<T, 3> uBoundary(param.getInletVelocity(0), 0.0, 0.0);
 
     if (param.lateralFreeSlip) {
         pcout << "Free-slip lateral boundaries on STL faces." << std::endl;
         lattice->periodicity().toggleAll(false);
         rhoBar->periodicity().toggleAll(false);
         j->periodicity().toggleAll(false);
 
         bc->setVelocityConditionOnBlockBoundaries(*lattice, param.inlet, boundary::dirichlet);
         setBoundaryVelocity(*lattice, param.inlet, uBoundary);
 
         bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral1, boundary::freeslip);
         bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral2, boundary::freeslip);
         bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral3, boundary::freeslip);
         bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral4, boundary::freeslip);
         setBoundaryVelocity(*lattice, param.lateral1, uBoundary);
         setBoundaryVelocity(*lattice, param.lateral2, uBoundary);
         setBoundaryVelocity(*lattice, param.lateral3, uBoundary);
         setBoundaryVelocity(*lattice, param.lateral4, uBoundary);
 
         Box3D globalDomain(lattice->getBoundingBox());
         std::vector<MultiBlock3D *> bcargs;
         bcargs.push_back(lattice);
         bcargs.push_back(rhoBar);
         bcargs.push_back(j);
         T outsideDensity = 1.0;
         int bcType = 1;
         integrateProcessingFunctional(
             new VirtualOutlet<T, DESCRIPTOR>(outsideDensity, globalDomain, bcType), param.outlet,
             bcargs, 2);
         setBoundaryVelocity(*lattice, param.outlet, uBoundary);
     } else {
         pcout << "Periodic lateral boundaries on STL faces." << std::endl;
         lattice->periodicity().toggleAll(true);
         rhoBar->periodicity().toggleAll(true);
         j->periodicity().toggleAll(true);
 
         lattice->periodicity().toggle(0, false);
         rhoBar->periodicity().toggle(0, false);
         j->periodicity().toggle(0, false);
 
         bc->addVelocityBoundary0N(param.inlet, *lattice);
         setBoundaryVelocity(*lattice, param.inlet, uBoundary);
 
         Box3D globalDomain(lattice->getBoundingBox());
         globalDomain.y0 -= 2;
         globalDomain.y1 += 2;
         globalDomain.z0 -= 2;
         globalDomain.z1 += 2;
         std::vector<MultiBlock3D *> bcargs;
         bcargs.push_back(lattice);
         bcargs.push_back(rhoBar);
         bcargs.push_back(j);
         T outsideDensity = 1.0;
         int bcType = 1;
         integrateProcessingFunctional(
             new VirtualOutlet<T, DESCRIPTOR>(outsideDensity, globalDomain, bcType), param.outlet,
             bcargs, 2);
         setBoundaryVelocity(*lattice, param.outlet, uBoundary);
     }
 }
 
 // Write VTK files for visualization
 void writeVTK(OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity> &bc, plint iT)
 {
     VtkImageOutput3D<T> vtkOut(createFileName("volume", iT, PADDING), param.dx);
     vtkOut.writeData<float>(
         *bc.computeVelocityNorm(param.boundingBox()), "velocityNorm", param.dx / param.dt);
     vtkOut.writeData<3, float>(
         *bc.computeVelocity(param.boundingBox()), "velocity", param.dx / param.dt);
     vtkOut.writeData<float>(
         *bc.computePressure(param.boundingBox()), "pressure",
         param.dx * param.dx / (param.dt * param.dt));
 }
 
 // Write PPM images on slices (adjusted for void center)
 void writePPM(OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity> &bc, plint iT)
 {
     // Approximate void center for slices (use domain center as fallback)
     plint cxLB = param.nx / 2;
     plint cyLB = param.ny / 2;
     plint czLB = param.nz / 2;
     Box3D xSlice(cxLB, cxLB, 0, param.ny - 1, 0, param.nz - 1);
     Box3D ySlice(0, param.nx - 1, cyLB, cyLB, 0, param.nz - 1);
     Box3D zSlice(0, param.nx - 1, 0, param.ny - 1, czLB, czLB);
 
     ImageWriter<T> writer("leeloo");
     writer.writeScaledPpm(
         createFileName("vnorm_xslice", iT, PADDING), *bc.computeVelocityNorm(xSlice));
     writer.writeScaledPpm(
         createFileName("vnorm_yslice", iT, PADDING), *bc.computeVelocityNorm(ySlice));
     writer.writeScaledPpm(
         createFileName("vnorm_zslice", iT, PADDING), *bc.computeVelocityNorm(zSlice));
 }
 
 
 void runProgram()
 {
     uint32_t seed = 1;
 
     /*
      * Read the domain geometry (STL with void as obstacle)
      */
     pcout << std::endl << "Reading STL data for the domain geometry." << std::endl;
     TriangleSet<T> triangleSet(param.geometry_fname, DBL);
     triangleSet.scale(1.0 / param.dx); // Convert to lattice units
     triangleSet.writeBinarySTL(outputDir + "domain_LB.stl");
 
     // Create boundary from STL
     plint xDirection = 0;
     plint borderWidth = 1;
     plint margin = 1;
     plint blockSize = 0;
     DEFscaledMesh<T> defMesh(triangleSet, 0, xDirection, margin, Dot3D(0, 0, 0));
     TriangleBoundary3D<T> boundary(defMesh);
 
     pcout << "tau = " << 1.0 / param.omega << std::endl;
     pcout << "dx = " << param.dx << std::endl;
     pcout << "dt = " << param.dt << std::endl;
     pcout << "Number of iterations in an integral time scale: " << (plint)(1.0 / param.dt)
           << std::endl;
 
     /*
      * Voxelize the domain (void as solid)
      */
     pcout << std::endl << "Voxelizing the domain." << std::endl;
     plint extendedEnvelopeWidth = 2;
     const int flowType = voxelFlag::outside; // Fluid outside void
     VoxelizedDomain3D<T> voxelizedDomain(
         boundary, flowType, param.boundingBox(), borderWidth, extendedEnvelopeWidth, blockSize);
     pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
 
     /*
      * Generate lattice and fields
      */
     pcout << "Generating the lattice, the rhoBar and j fields." << std::endl;
     MultiBlockLattice3D<T, DESCRIPTOR> *lattice =
         new MultiBlockLattice3D<T, DESCRIPTOR>(voxelizedDomain.getVoxelMatrix());
     if (param.useSmago) {
         defineDynamics(
             *lattice, lattice->getBoundingBox(),
             new SmagorinskyBGKdynamics<T, DESCRIPTOR>(param.omega, param.cSmago));
         pcout << "Using Smagorinsky BGK dynamics." << std::endl;
     } else {
         defineDynamics(
             *lattice, lattice->getBoundingBox(), new BGKdynamics<T, DESCRIPTOR>(param.omega));
         pcout << "Using BGK dynamics." << std::endl;
     }
     bool velIsJ = false;
     defineDynamics(
         *lattice, voxelizedDomain.getVoxelMatrix(), lattice->getBoundingBox(),
         new NoDynamics<T, DESCRIPTOR>(), voxelFlag::inside);
     lattice->toggleInternalStatistics(false);
 
     plint envelopeWidth = 1;
     MultiScalarField3D<T> *rhoBar =
         generateMultiScalarField<T>((MultiBlock3D &)*lattice, envelopeWidth).release();
     rhoBar->toggleInternalStatistics(false);
 
     MultiTensorField3D<T, 3> *j =
         generateMultiTensorField<T, 3>((MultiBlock3D &)*lattice, envelopeWidth).release();
     j->toggleInternalStatistics(false);
 
     std::vector<MultiBlock3D *> lattice_rho_bar_j_arg;
     lattice_rho_bar_j_arg.push_back(lattice);
     lattice_rho_bar_j_arg.push_back(rhoBar);
     lattice_rho_bar_j_arg.push_back(j);
     integrateProcessingFunctional(
         new ExternalRhoJcollideAndStream3D<T, DESCRIPTOR>(), lattice->getBoundingBox(),
         lattice_rho_bar_j_arg, 0);
     integrateProcessingFunctional(
         new BoxRhoBarJfunctional3D<T, DESCRIPTOR>(), lattice->getBoundingBox(),
         lattice_rho_bar_j_arg, 3);
 
     /*
      * Generate boundary conditions
      */
     pcout << "Generating boundary conditions." << std::endl;
 
     OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity> *boundaryCondition;
     BoundaryProfiles3D<T, Velocity> profiles;
     bool useAllDirections = true;
     OffLatticeModel3D<T, Velocity> *offLatticeModel = 0;
     if (param.freeSlipWall) {
         profiles.setWallProfile(new FreeSlipProfile3D<T>);
     } else {
         profiles.setWallProfile(new NoSlipProfile3D<T>);
     }
     offLatticeModel = new GuoOffLatticeModel3D<T, DESCRIPTOR>(
         new TriangleFlowShape3D<T, Array<T, 3> >(voxelizedDomain.getBoundary(), profiles), flowType,
         useAllDirections);
     offLatticeModel->setVelIsJ(velIsJ);
     boundaryCondition = new OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity>(
         offLatticeModel, voxelizedDomain, *lattice);
 
     boundaryCondition->insert();
 
     OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *outerBoundaryCondition =
         createLocalBoundaryCondition3D<T, DESCRIPTOR>();
     outerDomainBoundaries(lattice, rhoBar, j, outerBoundaryCondition);
 
     /*
      * Implement outlet sponge zone
      */
     if (param.numOutletSpongeCells > 0) {
         T bulkValue;
         Array<plint, 6> numSpongeCells;
 
         if (param.outletSpongeZoneType == 0) {
             pcout << "Generating an outlet viscosity sponge zone." << std::endl;
             bulkValue = param.omega;
         } else if (param.outletSpongeZoneType == 1) {
             pcout << "Generating an outlet Smagorinsky sponge zone." << std::endl;
             bulkValue = param.cSmago;
         } else {
             pcout << "Error: unknown type of sponge zone." << std::endl;
             exit(-1);
         }
 
         numSpongeCells[0] = 0;
         numSpongeCells[1] = param.numOutletSpongeCells;
         numSpongeCells[2] = 0;
         numSpongeCells[3] = 0;
         numSpongeCells[4] = 0;
         numSpongeCells[5] = 0;
 
         std::vector<MultiBlock3D *> args;
         args.push_back(lattice);
 
         if (param.outletSpongeZoneType == 0) {
             applyProcessingFunctional(
                 new ViscositySpongeZone3D<T, DESCRIPTOR>(
                     param.nx, param.ny, param.nz, bulkValue, numSpongeCells),
                 lattice->getBoundingBox(), args);
         } else {
             applyProcessingFunctional(
                 new SmagorinskySpongeZone3D<T, DESCRIPTOR>(
                     param.nx, param.ny, param.nz, bulkValue, param.targetSpongeCSmago,
                     numSpongeCells),
                 lattice->getBoundingBox(), args);
         }
     }
 
     /*
      * Set initial conditions
      */
     Array<T, 3> uBoundary(param.getInletVelocity(0), (T)0.0, (T)0.0);
     initializeAtEquilibrium(*lattice, lattice->getBoundingBox(), (T)1.0, uBoundary);
     applyProcessingFunctional(
         new BoxRhoBarJfunctional3D<T, DESCRIPTOR>(), lattice->getBoundingBox(),
         lattice_rho_bar_j_arg);
 
     /*
      * Particles (streamlines)
      */
     MultiParticleField3D<ParticleFieldT> *particles = 0;
     if (param.useParticles) {
         particles = new MultiParticleField3D<ParticleFieldT>(
             lattice->getMultiBlockManagement(),
             defaultMultiBlockPolicy3D().getCombinedStatistics());
 
         std::vector<MultiBlock3D *> particleArg;
         particleArg.push_back(particles);
 
         std::vector<MultiBlock3D *> particleFluidArg;
         particleFluidArg.push_back(particles);
         particleFluidArg.push_back(lattice);
 
         integrateProcessingFunctional(
             new AdvanceParticlesEveryWhereFunctional3D<T, DESCRIPTOR>(param.cutOffSpeedSqr),
             lattice->getBoundingBox(), particleArg, 0);
         integrateProcessingFunctional(
             new FluidToParticleCoupling3D<T, DESCRIPTOR>((T)param.particleTimeFactor),
             lattice->getBoundingBox(), particleFluidArg, 1);
 
         // Inject particles near inlet (adjusted for STL domain)
         Box3D injectionDomain(
             0, 0, param.ny / 4, 3 * param.ny / 4, param.nz / 4, 3 * param.nz / 4);
 
         Particle3D<T, DESCRIPTOR> *particleTemplate =
             new PointParticle3D<T, DESCRIPTOR>(0, Array<T, 3>(0., 0., 0.), Array<T, 3>(0., 0., 0.));
 
         std::vector<MultiBlock3D *> particleInjectionArg;
         particleInjectionArg.push_back(particles);
 
         integrateProcessingFunctional(
             new InjectRandomParticlesFunctionalPPRNG3D<T, DESCRIPTOR>(
                 particleTemplate, param.particleProbabilityPerCell, particles->getBoundingBox(),
                 &seed),
             injectionDomain, particleInjectionArg, 0);
 
         Box3D absorbtionDomain(param.outlet);
         integrateProcessingFunctional(
             new AbsorbParticlesFunctional3D<T, DESCRIPTOR>, absorbtionDomain, particleArg, 0);
 
         particles->executeInternalProcessors();
     }
 
     /*
      * Run simulation
      */
     plb_ofstream energyFile((outputDir + "average_energy.dat").c_str());
 
     pcout << std::endl;
     pcout << "Starting simulation." << std::endl;
     bool checkForErrors = true;
     for (plint i = 0; i < param.maxIter; ++i) {
         if (i <= param.initialIter) {
             Array<T, 3> uBoundary(param.getInletVelocity(i), 0.0, 0.0);
             setBoundaryVelocity(*lattice, param.inlet, uBoundary);
         }
 
         if (i % param.statIter == 0) {
             pcout << "At iteration " << i << ", t = " << i * param.dt << std::endl;
             if (i != 0) {
                 Array<T, 3> force(boundaryCondition->getForceOnObject());
                 T factor = util::sqr(util::sqr(param.dx)) / util::sqr(param.dt);
                 pcout << "Force on void (obstacle) over fluid density: F[x] = " << force[0] * factor
                       << ", F[y] = " << force[1] * factor << ", F[z] = " << force[2] * factor
                       << std::endl;
             }
             T avEnergy = boundaryCondition->computeAverageEnergy() * util::sqr(param.dx)
                          / util::sqr(param.dt);
             pcout << "Average kinetic energy over fluid density: E = " << avEnergy << std::endl;
             energyFile << i * param.dt << "  " << avEnergy << std::endl;
             pcout << std::endl;
         }
 
         if (i % param.vtkIter == 0) {
             pcout << "Writing VTK at time t = " << i * param.dt << endl;
             writeVTK(*boundaryCondition, i);
             if (param.useParticles) {
                 writeParticleVtk<T, DESCRIPTOR>(
                     *particles, createFileName(outputDir + "particles_", i, PADDING) + ".vtk",
                     param.dx, param.maxNumParticlesToWrite);
             }
         }
 
         if (i % param.imageIter == 0) {
             pcout << "Writing PPM image at time t = " << i * param.dt << endl;
             writePPM(*boundaryCondition, i);
         }
 
         lattice->executeInternalProcessors();
         lattice->incrementTime();
         if (param.useParticles && i % param.particleTimeFactor == 0) {
             particles->executeInternalProcessors();
         }
         seed++;
 
         if (checkForErrors) {
             abortIfErrorsOccurred();
             checkForErrors = false;
         }
     }
 
     energyFile.close();
     delete outerBoundaryCondition;
     delete boundaryCondition;
     if (param.useParticles) {
         delete particles;
     }
     delete j;
     delete rhoBar;
     delete lattice;
 }
 
 int main(int argc, char *argv[])
 {
     plbInit(&argc, &argv);
     global::directories().setOutputDir(outputDir);
 
     string xmlFileName;
     try {
         global::argv(1).read(xmlFileName);
     } catch (PlbIOException &exception) {
         pcout << "Wrong parameters; the syntax is: " << (std::string)global::argv(0)
               << " input-file.xml" << std::endl;
         return -1;
     }
 
     try {
         param = Param(xmlFileName);
     } catch (PlbIOException &exception) {
         pcout << exception.what() << std::endl;
         return -1;
     }
 
     try {
         runProgram();
     } catch (PlbIOException &exception) {
         pcout << exception.what() << std::endl;
         return -1;
     }
 }