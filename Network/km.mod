TITLE KA TOR
: K-A current for hippocampal interneurons from Lien et al (2002)


NEURON {
	SUFFIX km
	USEION k READ ek WRITE ik
	RANGE  gbar, ik, m1, h1, z1, z2, htau
	GLOBAL acthalf, actsteep, acttauavg, acttauoffset, acttauwdt, inactsteep, inactoffset, acttaupeak, inacttaupeak, inacttauxoffset, inacttauwdt, inacttauyoffset
}

PARAMETER {
	gbar = 0.0002   	(mho/cm2)	
								
	celsius
	ek		(mV)            : must be explicitly def. in hoc
	v 		(mV)
	z1
	z2
	q10=3
	acthalf = 50
	actsteep = 5
	acttauavg = 50
	acttauoffset = 57.5
	acttauwdt = 10
	acttaupeak = 10
	inactsteep = 15
	inactoffset = 30
	inacttaupeak = 450000
	inacttauxoffset = 57.5
	inacttauwdt = 30
	inacttauyoffset = 9
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	minf 		mtau (ms)
	hinf	 	htau (ms)
}
 

STATE { m1 h1 }

BREAKPOINT {
        SOLVE states METHOD cnexp
	ik = gbar*m1*h1*(v - ek)
	z1 = m1*h1
	z2 = m1*h1*(v - ek)
} 

INITIAL {
	trates(v)
	m1=minf  
	h1=hinf	
}

DERIVATIVE states {   
        trates(v)      
        m1' = (minf-m1)/mtau
		h1' = (hinf-h1)/htau
}

PROCEDURE trates(v) {  

	
	LOCAL qt
        qt=q10^((celsius-23)/10)
		:Costum functions for better fitting
		
        minf = (1/(1 + exp(-(v+acthalf)/actsteep)))^4
	mtau=acttauavg+acttaupeak*exp(-0.5*((v+acttauoffset)/acttauwdt)^2)
		hinf = 1/(1 + exp((v+inactoffset)/inactsteep))
		htau = inacttauyoffset+inacttaupeak*exp(-0.5*((v+inacttauxoffset)/inacttauwdt)^2)
	


}