! --------------------------------------------------------------------------
! File: pp13TeV_ttbar.cmd
! This file contains commands to be read in for a Pythia8 run. 
! Lines not beginning with a letter or digit are comments.
! Names are case-insensitive  -  but spellings-sensitive!
! The changes here are illustrative, not always physics-motivated.
! --------------------------------------------------------------------------
! 1) Settings that will be used in a main program.
Main:numberOfEvents =  200000      ! number of events to generate/read
Main:timesAllowErrors = 10000      ! abort run after this many flawed events
! --------------------------------------------------------------------------
! 2) Settings related to output in init(), next() and stat().
Init:showChangedSettings = on      ! list changed settings
Init:showAllSettings = off         ! list all settings
Init:showChangedParticleData = on  ! list changed particle data
Init:showAllParticleData = off     ! list all particle data

Next:numberCount = 10000           ! print message every n events
Next:numberShowLHA = 0             ! print LHA information n times
Next:numberShowInfo = 0            ! print event information n times
Next:numberShowProcess = 0         ! print process record n times
Next:numberShowEvent = 0           ! print event record n times

Stat:showPartonLevel = off         ! additional statistics on MPI
! --------------------------------------------------------------------------
! 3) Beam parameter settings. 
Beams:idA = 2212                   ! first beam,  p = 2212, pbar = -2212
Beams:idB = 2212                   ! second beam, p = 2212, pbar = -2212
Beams:eCM = 13000.                 ! CM energy of collision
! --------------------------------------------------------------------------
! 4) Process
Top:all     = on		   ! all production modes
25:m0       = 125.0		   ! (GeV) Higgs boson mass

PhaseSpace:pTHatMin = 2.	   ! minimum pT of hard process
! --------------------------------------------------------------------------
! 5) Parton level control.
PartonLevel:ISR = on               ! Initial state radiation
PartonLevel:FSR = on               ! Final state radiation
PartonLevel:MPI = on		   ! Multiple interactions
! --------------------------------------------------------------------------
Random::setSeed = on
Random:seed     = 4494
