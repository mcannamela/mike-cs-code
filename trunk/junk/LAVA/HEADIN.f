C=============================================================HEADING
C
C     LAVA \LA-VA\ N. 1: FLUID ROCK THAT ISSUES FROM A VOLCANO OR
C                     FROM A FISSURE IN THE EARTH'S SURFACE, ALSO
C                     SUCH ROCK SOLIDIFIED.
C                     2: A COMPUTER PROGRAM FOR TWO OR THREE
C                     DIMENSIONAL FLUID FLOW WITH CHEMICAL REACTIONS
C                     AND DISPERSED PARTICLES.
C
C     BY    C. H. CHANG  AND J. D. RAMSHAW
C           IDAHO NATIONAL ENGINEERING LABORATORY
C
C           JAN. 20, 1997
C
C     THE DEVELOPMENT OF THIS PROGRAM HAS BEEN SUPPORTED BY THE U. S.
C     DEPARTMENT OF ENERGY, OFFICE OF BASIC ENERGY SCIENCES (BES)
C     PROGRAM.
C
C=============================================================HEADING
C
C     COPYRIGHT
C
C=============================================================HEADING
C     THIS MATERIAL RESULTED FROM WORK DEVELOPED UNDER COVERNMENT
C     CONTRACT NO. DE-AC07-94ID13223 AT THE IDAHO NATIONAL
C     ENGINEERING LABORATORY AND IS SUBJECT TO A BROAD GOVERNMENT
C     LICENSE.  NEITHER THE UNITED STATES NOR THE UNITED STATES
C     DEPARTMENT OF ENERGY, NOR ANY OF THEIR EMPLOYEES MAKE ANY
C     WARRANTY, EXPRESS OR IMPLIED, INCLUDING THE WARRANTIES OF
C     MERCHANTABILITY OR FITNESS FOR A PRTICULAR PURPOSE, OR ASSUMES
C     ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
C     COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS,
C     PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD
C     NOT INFRINGE PRIVATELY-OWNED RIGHTS.
C=============================================================HEADING
C
C     FACE AND INDEX NOTATION
C
C                                         _____________________ NXYZT
C                                        /|                   /|
C                                      /                    /  |
C                                    /    |  TOP          /    |
C                                  /             FRONT  /      |
C                                /        |           /        |
C                              /____________________/          |
C        Z                     |   LEFT   |_ _ _ _ |_ _RIGHT_ _|
C        |        Y            |                   |          / NXYT
C        |      /              |        /          |        /
C        |    /                |      DERRIERE     |      /
C        |  /                  |    /       BOTTOM |    /
C        |/                    |                   |  /
C        ----------> X         |/__________________|/
C                             1                     NXT
C
C=============================================================HEADING
C     THIS PROGRAM IS WRITTEN IN DOUBLE PRECISION.
C     FOR CRAY, COMPILER OPTION -DP (LOWER CASES) SHOULD BE USED TO
C     IGNORE DOUBLE PRECISION.
C=============================================================HEADING
C     SHOPPING LIST FOR LAVA
C-------------------------------------------------------------HEADING
C     radiative transport
C-------------------------------------------------------------HEADING
C     Particle size distribution -- modify so that sieve data with
C     particle size distribution interval can be input directly.
C-------------------------------------------------------------HEADING
C     multiple particle type, multiple injection
C     evaporation model
C-------------------------------------------------------------HEADING
C     multicompoent particle (condensed phase logic) -- modify both
C      INJECT and PSTATE
C-------------------------------------------------------------HEADING
C     CHMIMP -- include three-body reactions (chaffron efficiency)
C     Chemistry source/sink due to particles
C     Kinetic reaction for heterogeneous reactions
C      (Currently, only equilibrium in CHMIMP is allowed.)
C     To handle reactions with
C      Ke (T,Te) = kf(Te)/kb(T)  (e.g. ArH+ + e- <--> H* + Ar )
C     Dynamically switch forward and backward reaction based on the
C      PP, RP, and EQC.  Do not do this in the middle of iteration.
C-------------------------------------------------------------HEADING
C     Combine input variables for both implicit and explicit chemistry
C      routines.  Invent input flag to tell whether reaction is fast
C      or slow.
C-------------------------------------------------------------HEADING
C     Implement and test nucleation and condensation in stochastic
C      particle model.  Compute mass and energy exchange using
C      chemistry routines.
C-------------------------------------------------------------HEADING
C     Diffusion -- include thermal and pressure diffusion
C               -- implement friction-weighted approximation
C               -- implement Lennard-Jones parameters
C-------------------------------------------------------------HEADING
C     Check whether volumetric radiation energy loss data includes
C     resonance radiation.  If so, that term needs to be reduced by
C     80 % or so.  See Jim Menart's work
C-------------------------------------------------------------HEADING
C     Finalize wall function and SETR for curved geometries
C      for ellipses, straight walls, spheres, cylinders, etc.
C     Generalize SCORELS obstacle logic for cylindrical coord.
C-------------------------------------------------------------HEADING
C     Terms for Time-dependent R factors
C-------------------------------------------------------------HEADING
C     Turbulent diffusion of condensed species (currently it is
C      assumed that eddy motions (consequently diffusion coefficients
C      are the same for both gas and condensed phases).
C     Develop APACHE like formula for mixture of gas and condensed
C      phases, instead of multiplying by IGAS.
C=============================================================HEADING
C     Multitemperature multicomponent diffusion with superimposed EM
C     field
C-------------------------------------------------------------HEADING
C     Vibratrional non-equilibrium and its effects on charge exchnage
C-------------------------------------------------------------HEADING
C     Minimization technique for kinetic and equilibrium chemistry
C      Combine Gibbs free energy and kinetic rate and minimize for
C      general purpose chemistry routine - sell
C=============================================================HEADING
