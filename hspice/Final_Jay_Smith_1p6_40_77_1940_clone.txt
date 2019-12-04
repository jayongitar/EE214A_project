* Design Project EE 114/214A - 2019
* Team Member 1 Name:  Jay Smith
* Team Member 2 Name:
* Please fill in the specification achieved by your circuit 
* before you submit the netlist.
**************************************************************
* SUNET IDs of team members = 06407853
* The specifications that this script achieves are: 
* Power     = 1.590 mW
* Gain      = 40.074 kohms
* Bandwidth = 76.926 MHz
* FOM       = 1939
* Bonus
* Spot Noise @ 1 kHz = 
* FOM =  
***************************************************************
.subckt tia iina iinb vouta voutb vdd vss

** Including the model file
.include /afs/ir.stanford.edu/class/ee114/hspice/ee114_hspice.sp

*** Your Transimpedance Amplifier
** parameters 
.param ml1_w = 5.2u    ml1_l = 2u
.param m1_w =  6.4u    m1_l  = 1u
.param mb1_w = 4.2u    mb1_l = 2u

.param ru = 57k	  rd = 73k

.param ml2_w = 3.2u    ml2_l = 2u
.param m2_w =  7.0u    m2_l = 1u
.param mb2_w = 2.8u    mb2_l = 2u

.param m3_w  = 36.6u   m3_l = 1u
.param mb3_w = 9.6u    mb3_l = 2u

.param iref  = 20u
.param mi1_w = 3.2u   mi1_l = 2.0u
.param mi2_w = 3.2u   mi2_l = 2.0u
.param mi3_w = 4.4u   mi3_l = 2.0u
**	d	g	s	b	n/pmos114	w	l
*** A Side ***

* CG stage 
ml1a	cg_outa	  vb_pmos	vdd	    vdd	 pmos114	 w='ml1_w'	l='ml1_l'
m1a	    cg_outa	  0	        iina	vss	 nmos114	 w='m1_w'	l='m1_l'
mb1a	iina	  vb_nmos	vss	    vss	 nmos114	 w='mb1_w'	l='mb1_l'
rua	    vdd	      cg_outa	'ru'
rda	    cg_outa	  vss	    'rd'

* CS stage
ml2a	vdd	      vdd	    cs_outa	vss	 nmos114	w='ml2_w'	l='ml2_l'
m2a	    cs_outa	  cg_outa	d_gnd	vss	 nmos114	w='m2_w'	l='m2_l'
mb2a	d_gnd	  vb_nmos	vss	    vss	 nmos114	w='mb2_w'	l='mb2_l'

* CD stage
m3a	    vdd	      cs_outa   vouta	vss	 nmos114	w='m3_w'	l='m3_l'
mb3a	vouta	  vb_nmos	vss	    vss	 nmos114    w='mb3_w'	l='mb3_l'



*** B Side ***
* CG stage 
ml1b	cg_outb	vb_pmos	vdd	vdd	pmos114	w='ml1_w'	l='ml1_l'
m1b	cg_outb	0	iinb	vss	nmos114	w='m1_w'	l='m1_l'
mb1b	iinb	vb_nmos	vss	vss	nmos114	w='mb1_w'	l='mb1_l'
rub	vdd	cg_outb	'ru'
rdb	cg_outb	vss	'rd'

* CS stage
ml2b	vdd	vdd	cs_outb	vss	nmos114	w='ml2_w'	l='ml2_l'
m2b	cs_outb	cg_outb	d_gnd	vss	nmos114	w='m2_w'	l='m2_l'
mb2b	d_gnd	vb_nmos	vss	vss	nmos114	w='mb2_w'	l='mb2_l'

* CD stage
m3b	vdd	cs_outb voutb	vss	nmos114	w='m3_w'	l='m3_l'
mb3b	voutb	vb_nmos	vss	vss	nmos114 w='mb3_w'	l='mb3_l'

*** Current Bias ***
iref    vdd      vb_nmos  'iref'
mi1     vb_nmos  vb_nmos  vss     vss     nmos114 w='mi1_w'       l='mi1_l'
mi2     vb_pmos  vb_nmos  vss     vss     nmos114 w='mi2_w'       l='mi2_l'
mi3     vb_pmos  vb_pmos  vdd     vdd     pmos114 w='mi3_w'       l='mi3_l'


.ends