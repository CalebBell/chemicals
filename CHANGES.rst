




BREAKING CHANGES FROM THERMO
* Many functions and classes are not included here; they will remain in `thermo`
* The optional "available methods" behavior for data retrieval functions has been refactored to their own specialized function. For example, to get available critical temperature methods for a given CASRN, use `Tc_methods(CASRN)` instead of `Tc(CASRN, AvailableMethods=True)`.

BREAKING CHANGES IN THINGS ADDED SINCE DEVELOPMENT OF CHEMICALS STARTED
*