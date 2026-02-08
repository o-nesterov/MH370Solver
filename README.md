# MH370Solver
MH370 path modelling software

The flight path modelling software developed to assist in finding the terminal location of the Malaysian MH370, Boeing-777, 9M-MRO, which went missing on 8 March 2014.

The approach is based on the automated minimization of the burst timing offset and burst frequency offset differences between the model and Inmarsat measurements by varying user-selected flight path parameters, such as the initial location, altitude, IAS/Mach, duration of legs, etc. A flight path is comprised of the set of maneuvers, presently legs and turns, defined by a user. The path is obtained by numerical integration of the ODE systems using the 5(4) RK method with an adaptive integration time step. The fuel flow model is included. The meteorological forcing is optionally sourced either from GDAS1 or ERA5, with several interpolation methods incorporated (e.g., linear, cubic, splines). Several heading modes are supported: true track, true heading, magnetic track, magnetic heading, and gyroscopic heading. Two speed modes are supported: IAS and Mach. Supported vertical modes include level flight, V/S and FPA level change.

This project is not owned by, affiliated with, or sponsored by any governmental, commercial or private organization.
