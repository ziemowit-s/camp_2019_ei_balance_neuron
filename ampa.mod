TITLE minimal model of AMPA receptors

COMMENT
-----------------------------------------------------------------------------

	Minimal kinetic model for glutamate AMPA receptors
	==================================================

  Model of Destexhe, Mainen & Sejnowski, 1994:

	(closed) + T <-> (open)

  The simplest kinetics are considered for the binding of transmitter (T)
  to open postsynaptic receptors.   The corresponding equations are in
  similar form as the Hodgkin-Huxley model:

	dr/dt = alpha * [T] * (1-r) - beta * r

	I = gmax * [open] * (V-Erev)

  where [T] is the transmitter concentration and r is the fraction of 
  receptors in the open form.

-----------------------------------------------------------------------------

  Based on voltage-clamp recordings of AMPA receptor-mediated currents in rat
  hippocampal slices (Xiang et al., J. Neurophysiol. 71: 2552-2556, 1994), this
  model was fit directly to experimental recordings in order to obtain the
  optimal values for the parameters (see Destexhe, Mainen and Sejnowski, 1996).

-----------------------------------------------------------------------------

  This mod file includes a mechanism to describe the time course of transmitter
  on the receptors.  The time course is approximated here as a brief pulse
  triggered when the presynaptic compartment produces an action potential.
  The trigger is sensed by the NET_RECEIVE function attached at the
  last, which sets the value of C to Cmax once the trigger is received.

-----------------------------------------------------------------------------

  See details in:

  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for
  computing synaptic conductances based on a kinetic model of receptor binding
  Neural Computation 6: 10-14, 1994.  

  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of 
  synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition; 
  edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1998, pp. 1-25.

-----------------------------------------------------------------------------
ENDCOMMENT



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON 
{
	POINT_PROCESS AMPA
	RANGE C, g, gmax, lastrelease, TRise, tau
	NONSPECIFIC_CURRENT i
	RANGE Cmax, Cdur, Alpha, Beta, Erev, Deadtime
}

UNITS 
{
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER 
{
	TRise (ms)
	tau   (ms)
	Cmax	= 1	(mM)		: max transmitter concentration
	Erev	= 0	(mV)		: reversal potential
	Deadtime = 1	(ms)		: mimimum time between release events
	gmax		(umho)		: maximum conductance
}


ASSIGNED 
{
	Alpha	(/ms mM)	: forward (binding) rate
	Beta	(/ms)		: backward (unbinding) rate
	Cdur	(ms)		: transmitter duration (rising phase)
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(umho)		: conductance
	C		(mM)		: transmitter concentration
	lastrelease	(ms)		: time of last spike
}

STATE
{
	R				: fraction of open channels
}

INITIAL 
{
	R = 0
	C = 0
	lastrelease = -1000
	Cdur=TRise
	Beta=1/tau
	Alpha=1/Cdur - Beta
}

BREAKPOINT 
{
	SOLVE states METHOD euler
	g = (gmax * R * (Alpha+Beta)) / (Alpha*(1-1/exp(1)))
	i = g*(v - Erev)
}

DERIVATIVE states
{
	evaluateC() 	: Find out value of C
	R'=Alpha * C * (1-R) - Beta * R
}

PROCEDURE evaluateC()
{
	LOCAL q
	q = ((t - lastrelease) - Cdur)		: time since last release ended
	if (q >= 0 && q <= Deadtime && C == Cmax) {	: in dead time after release
		C = 0.
	}
}

NET_RECEIVE (weight (umho)) 
{ 
	LOCAL q
	q = ((t - lastrelease) - Cdur)		: time since last release ended

: Spike has arrived, ready for another release?

	if (q > Deadtime) {
		C = Cmax			: start new release
		lastrelease = t
	} 
}

