## SDPT3 version 4.0: MATLAB/Octave software for semidefinite-quadratic-linear programming
### [Kim-Chuan Toh](https://blog.nus.edu.sg/mattohkc/), [Michael J. Todd](https://people.orie.cornell.edu/miketodd/todd.html), and [Reha H. Tütüncü](https://www.math.cmu.edu/~reha/)

SDPT3 is a Matlab package for solving convex optimization problems involving linear equations and inequalities, second-order cone constraints, and semidefinite constraints (linear matrix inequalities).

This is an *unofficial* repository for SDPT3. The [official SDPT3 site](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/) and [another GitHub repository](https://github.com/Kim-ChuanToh/SDPT3) is administered by [Kim-Chuan Toh](https://blog.nus.edu.sg/mattohkc/), co-author of SDPT3 and professor of mathematics at the National University of Singapore.

This repo is administered by [Michael Grant](http://cvxr.com/bio), the developer of [CVX](http://cvxr.com/cvx), a modeling framework for convex optimization that uses SDPT3 as a solver. *It is not our intent for this repo to become an independent fork of SDPT3*. To that end:

   + We will periodically monitor the [official site](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/) for updates, and incorporate any we find into this repo. 
   + Any improvements that we make will be submitted to the authors for inclusion in the official release.

Having said that, this distribution *does* differ from the release currently posted on the main site in several ways:

   + *Preliminary* support for Octave 3.8.0 and later has been added. While we have found that it works well,
     you should expect an occasional error or incompatibility simple due to a lack of testing.
   + The installer has been replaced. The new installer will build the MEX files only if they are not already
     present (or unless a "-rebuild" flag is given). It will also set the MATLAB paths, if necessary.
     *At this point, Octave users will need to build their own MEX files.*
   + The ".dll" files have been removed. Those using very old versions of MATLAB will need to rebuild them.
     The maintainer no longer keeps a version of MATLAB that old, so they can no longer be verified.
   + The code has been reformatted and cleaned up to eliminate all warnings from MATLAB's Code Analyzer.
   + Old files that are no longer in use have been removed. (That's what a repo is for!)

While you are welcome to submit bug reports on the [GitHub issue page](https://github.com/sqlp/sdpt3/issues) for this repo, we cannot guarantee that they will be addressed in a timely fashion.

Citation:

   + K.C. Toh, M.J. Todd, and R.H. Tutuncu, SDPT3 — a Matlab software package for semidefinite programming, Optimization Methods and Software, 11 (1999), pp. 545–581.
   + R.H Tutuncu, K.C. Toh, and M.J. Todd, Solving semidefinite-quadratic-linear programs using SDPT3, Mathematical Programming Ser. B, 95 (2003), pp. 189–217.

This version of SDPT3 is distributed under the GNU General Public License 2.0. For more details, please see the included files [Copyright](https://github.com/sqlp/sdpt3/blob/master/Copyright) and [GNU\_General\_Public\_License\_v2](https://github.com/sqlp/sdpt3/blob/master/GNU_General_Public_License_v2).

Michael C. Grant   
[CVX Research, Inc.](http://cvxr.com)
