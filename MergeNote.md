# Merge note of macana svn vs cvs

- Analysis/AnalParams.cpp
	- what does macanaRadom do?
		- for now I'm keeping it
	- remove simparam
	- add icamode (int)
	- add mask (int)
	- remove altKernel
	- remove simulate, subtract, and subPath
	- ! NcFile constructor
- Analysis/SimParams.cpp
	- xMap?
	- add atmFreq resample/atmSeed
- Clean/AzElTemplateCalculator.cpp
	- add tCoeffs
	- fix npars (minus 1)
	- add downWeight (What does this do)?
	- add CUBIC mode
	- why suspicious gsl_vector_scale(tmp, 0.9)?
	- add method calculateTemplate2
		- this is actually used by the azelResidual, why?
		- what is the difference in between this one and calcualteTemplate?
		- some suspicious stuff going on
			- why check on number of masked detectors: 3/4 is the magic number?
	- orphan method decorrelateCoeffs
		- why multiple by 0.9?
	- add CUBIC to setMode (via calculateMode)
		- replace AVERAGE with CONST?
	- change overrideMode
		- !!this implementation does not use the emum AzElFitType but uses int
	- add updateTemplate
- Clean/Clean.cpp
	- add nodc=false (decorrelate?)
- Clean/CleanBspline.cpp
	- add macros USEMPFIT=0 and LINFIT=1
	- set bright using ap->mask in CleanBspline
		- bright in degrees
		- ap->mask in arcsec
		- was 0.5arcmin for LMT, was 2arcmin for default 
	- add cleanKernel=false
	- add refBolo in removeBadBolos
	- change from mean to linfit_flags
	- change 3.0 to 3.5 for detPow/mPow for badBolo
	- add conditional remove badbolo ap->getOrder != 100 ap->getAzElMap !=1
	- return immediately when ap->getOrder == 100 (100 means atomosphere clean turned off by user)
	- cleanKernel controls wether subtract atmTempalte from detector.hKernel 
	- add clean-up resources at the end of clean
		- why do removeBadBolos again after destroyBaseMatrix?
	- add check on valid resample value
	- change cleanScans loop
		- azElRaPhys to hAzPhys
		- azElDecPhys to hElPhys
		- bypass hSampleFlag to all true ??
		- dist now in radians, and convert to degree
	- fixFalgs condition
		- cottingham use refbolo (-1) instead of 0
	- remove STRIPE_SPLINE and STRIPE_SCAN from StripeMethod
	- cottingham
		- make lin_fit optional. use mean as alternative
		- remove engergy and meanTod
		- add cleanup and return NULL on nan detected in tsi vector; instead of exit
		- remove pointDist check on bright
	-  remove removeScanPattern
	- change detSens to unit of Jy
		- was doing 1/detector.getSensitivity, why???
	- add flagVectorDist for flags of the bright mask
	- change removeLargeScaleResiduals loop
		- azElRaPhys to hAzPhys
		- azElDecPhys to hElPhys
		- dist with hRA and hDec
	- check detSens isnan. if so exit due to no good data on Az El template
	- update flagVector to detector.hSample
	- add createBaseMatrix overrivde with azOff and elOff
		- Need find out the different of this ??
	- remove removeCorrelation and removeCorrelation2
	- rewrote pca from stripPCA
	- add check of sample value in subtractTemplate
	- remove cleanStripeFFT2
	- change cleanStripeFFT3 to cleanStripFFT
	- remove splineResidual
	- use calculateTemplate2 in azelResidualLMT
	- use ap->cleanStripe to override clean mode
	- set corrCoefs to 1/sens ** 2
	- remove AzelResidual
- Clean/CleanHigh.cpp
	-  loops with  nrow x nrow. need check?
- Clean/CleanPCA.cpp
	- implement ICAmode when neigToCut >0
- Clean/CleanSelector.cpp
	- change order>3 from CleanHigh to CleanIterativeMean
- MapMaking/Coaddition.cpp
	- write noise to coadd map
	- add metadata nDataFiles and totalIntTime
	- add ap keywords icamode, bsplineMethod, and bsplineMaks
- MapMaking/CompletenessSim.cpp
	- use GslRandom instead of macanaRandom
	- is ap->macanaRandom obsolete? need check
- MapMaking/Map.cpp
	- !!diff failed. for now using theirs. will revisit later
- MapMaking/Observation.cpp
	- mapGenerationPrep
		- !!check for consistency
- Mapmaking/WienerFilter
	- done when maxtratio <1e-10
- Mapmaking/Array.cpp
	- fakeAtmData atmFreq is now set to 0.5  (was 0.05)
	- add refBoloIndex logic to handle Muscat
- Observatory/Detector.cpp
	- change fix hSampleFlags loop index j from 0 to 1 according to GW's note
	- Keep setFCF function
- Observatory/Telescope.cpp
	- add check on turning signal to prevent large map size
- Sky/Novas
	- add thread safety code 
- Utilities/BinomialStats.cpp
	- fix a typo of assigning np in calcIntervals
- Utilities/vector_utilites.cpp
	- !!diff failed. for now using theirs. will revisit later 
-  Utilities/vector_utilites.h
	- !!diff failed. for now using theirs. will revisit later 



















