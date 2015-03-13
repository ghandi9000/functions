## Functions used in bole volume and height prediction

##################################################################
######						                                ######
######		Formulas for bole volume calculation			######
######						                                ######
##################################################################
### kozak taper equation, returns diameter, d, at specified height, h
### funny volumes below 0.9 DBH, so small trees excluded
kozak2004 <- function(D, H, species, h) {
	if(species == "ABBA") {
		a0 = 0.911
		a1 = 1.026
		a2 = -0.005
		b1 = 0.368
		b2 = -0.645
		b3 = 0.502
		b4 = 1.780
		b5 = 0.096
		b6 = -0.487
	}
	if(species == "PIRU") {
		a0 = 0.940
		a1 = 0.998
		a2 = 0.010
		b1 = 0.508
		b2 = -0.636
		b3 = 0.355
		b4 = 1.687
		b5 = 0.078
		b6 = -0.242
	}
	z <- h/H
	p <- 1.3/H
	Q <- 1-z^(1/3)
	X <- (1-(h/H)^(1/3))/(1-p^(1/3))
	d <- a0*(D^a1)*(H^a2)*X^(b1*(z^4)+b2*(1/exp(D/H))+b3*(X^0.1)+b4*(1/D)+b5*(H^Q)+b6*X)
	return( d );
}

### Smalian formula, D1 and D2 are diameters at each end of log, L is length
### returns volume in cubic m
smalian <- function(D1, D2, L) {
	D1 <- 0.00007854*(D1^2)
	D2 <- 0.00007854*(D2^2)
	volume = ((D1+D2)/2)*L
	return( volume )
}

### Honer eq 1965
honer <- function(D, H, species) {
	D <- 0.393700787*D
	H <- 3.2808399*H
	if(species == "ABBA") {
		a = 2.139
		b = 301.634
	}
	if(species == "PIRU") {
		a = 0.691
		b = 363.676
	}
	V = 0.0283168466*D^2/(a+(b/H))
	return( V )
}

### this version of the honer equation takes dbh in cm and height in m
### thus returning cubic meter
honer2 <- function(dbh, height, species) {
	if(species == "BECO") { a0 = 2.222; a1 = 91.554; a2 = 0.0043222 }
	if(species == "BEAL") { a0 = 1.449; a1 = 105.081; a2 = 0.004320 }
	if(species == "ACSA" | species == "ACSP" | species == "ACPE") {
		a0 =1.046; a1 = 117.035; a2 = 0.004334 }
	if(species == "FAGR") {
		d1.3 <- dbh/(1-0.04365*0.145);
		dbh <- d1.3;
		a0=0.959; a1 = 102.056; a2 = 0.004334 }
	if(species == "PRPE" | species == "SOAM") { a0 = 0.033; a1 = 119.889; a2 = 0.004334 }
	V = (a2*dbh^2)/(a0+(a1/height))
	return( V )
}

###################### Bole Volume Formula For Kozak BV ##########################
bolevol <- function(dbh, ht, species) {
	increment <- ht/10
	heights <- seq(from=0,to=ht, by = increment)
	diams <- kozak2004(dbh,ht,species=species,heights)
	volume = 0
	for(j in 1:10) { volume=volume+smalian(diams[[j]],diams[[j+1]],increment) }
	return( volume )
}

###################### Clark BV Equation #########################################
### Requires trees to be 17.3 feet tall and 5 in diameter
### input units are dbh in cm, height in m, output is in m^3

clark <- function(dbh, height, species) {
	##### Species Parameters ########################
	if(species=="PRPE" | species=="SOAM") { a1 = 0.92487; b1 = -0.89867;  r = 37.12714; c = 0.48776; e = 1.50579;
		p = 6.18866; b = 1.64261; a = 0.55071 }
###	if(species=="ABBA" | species=="PIRU") { a1 = 0.92487; b1 = -0.89867;  r = 31.66250; c = 0.57402; e = 110.96;
###		p = 8.573; b = 2.36238; a = 0.68464 } ### using the coefficients for loblolly pine to get comparisons for ABBA and PIRU

	if(species=="ACPE" | species=="ACSA" | species=="ACSP") { a1 = 0.93991; b1 = -1.62226; r = 22.00135; c = 0.45472;
		e = 166.1; p = 7.31546; b = 1.17064; a = 0.27213 }
	if(species=="FAGR") { a1 = 0.91141; b1 = -0.6673; r = 44.36826; c = 1.22158; e = 79.44636;
		p = 6.36236; b = 1.11382; a = 0.14312 }
	if(species=="BEAL" | species=="BECO" | species=="BEPA") { a1 = 0.85516; b1 = -0.00134; r = 49.41385; c = 1.01241;
		e = -91.82769; p = 11.23179; b = 1.19704; a = 0.23928 }
	###if(species=="SOAM") { a1 = 0.91531; b1 = -0.96788; r = 0.64600; c = 0.49680; e = 127.87;
	###	p = 7.52915; b = 1.49287; a = 0.47222 }

	### Basic Symbols ###############################
	### V <- stem volumen between L and U in cubic ft.
	### L <- lower height of interest
	### U <- upper height of interest
	### D <- diameter in inches at 4.5 ft
	### H <- total ht of tree in FT
	### F <- diameter at 17.3 ft above ground in inches, DOB17 = D(a+b(17.3/H)^2) where a,b are sp. coefs
	D <- dbh*0.393700787 ### convert dbh in cm to dbh in in.
	L <- 0
	U <- height*3.2808399 ### convert height in meters to height in feet
	H <- U
	F <- D*(a1+b1*(17.3/H)^2)
	### Combined Variables #########################
	L1 <- max(L,0)
	U1 <- min(U,4.5)
	L2 <- max(L,4.5)
	U2 <- min(U,17.3)
	L3 <- max(L,17.3)
	U3 <- min(U,H)
	G <- (1-4.5/H)^r
	W <- (c + e/D^3)/(1-G)
	X <- (1-4.5/H)^p
	Y <- (1-17.3/H)^p
	Z <- (D^2-F^2)/(X-Y)
	T <- (D^2-Z*X)

	### Indicator Variables ########################
	I1 <- if(L < 4.5) I1 <- 1 else(I1 <- 0)
	I2 <- if(L < 17.3) I2 <- 1 else(I2 <- 0)
	I3 <- if(U > 4.5) I3 <- 1 else(I3 <- 0)
	I4 <- if(U > 17.3) I4 <- 1 else(I4 <- 0)
	I5 <- if((L3-17.3)<a*(H-17.3)) I5 <- 1 else(I5 <- 0)
	I6 <- if((U3-17.3)<a*(H-17.3)) I6 <- 1 else(I6 <- 0)

	### main stem volume calculation ################
	V<- 0.005454154*(I1*D^2*((1-G*W)*(U1-L1)+W*((1-L1/H)^r*(H-L1) -
	(1-U1/H)^r*(H-U1))/(r+1))
	+ I2*I3*(T*(U2-L2)+Z*((1-L2/H)^p*(H-L2) -
	(1-U2/H)^p*(H-U2))/(p+1))
	+ I4*F^2*(b*(U3-L3)-b*((U3-17.3)^2-(L3-17.3)^2)/(H-17.3) +
	(b/3)*((U3-17.3)^3-(L3-17.3)^3)/(H-17.3)^2 +
	I5*(1/3)*((1-b)/a^2)*(a*(H-17.3)-(L3-17.3))^3/(H-17.3)^2 -
	I6*(1/3)*((1-b)/a^2)*(a*(H-17.3)-(U3-17.3))^3/(H-17.3)^2))
	return( V*0.0283168466 )
}

## This is awful, but its the function used to create subsets of data used during
## height calculation in htsFINAL.R
make.working.dataset <- function(spp) {
	trees.from.86 <- which ( 	pp$SPECIES == spp & !is.na(pp$HT86) & !is.na(pp$DBH86) & pp$DBH86 != 0 )
	trees.from.87 <- which ( 	pp$SPECIES == spp & !is.na(pp$HT87) & !is.na(pp$DBH87) & pp$DBH87 != 0 )
	trees.from.87 <- trees.from.87[is.na(match(trees.from.87,trees.from.86))]
	trees.from.98 <- which (	pp$SPECIES == spp & !is.na(pp$HT98) & !is.na(pp$DBH98) & pp$DBH98 != 0 )
	trees.from.98 <- trees.from.98[is.na(match(trees.from.98,trees.from.86))]
	trees.from.98 <- trees.from.98[is.na(match(trees.from.98,trees.from.87))]
	trees.from.10 <- which (	pp$SPECIES == spp & !is.na(pp$HT10) & !is.na(pp$DBH10) & pp$DBH10 != 0 )
	trees.from.10 <- trees.from.10[is.na(match(trees.from.10,trees.from.86))]
	trees.from.10 <- trees.from.10[is.na(match(trees.from.10,trees.from.87))]
	trees.from.10 <- trees.from.10[is.na(match(trees.from.10,trees.from.98))]
	dbh86 <- pp[trees.from.86,]$DBH86
	ht86 <- pp[trees.from.86,]$HT86
	elev86 <- pp[trees.from.86,]$ELEV
	plots86 <- pp[trees.from.86,]$PLOT
	canht86 <- pp[trees.from.86,]$CANHT86
	year86 <- rep(86,length(ht86))
	data86 <- data.frame(plots86,dbh86,ht86,canht86,elev86,year86)
	names(data86) <- c("PLOT","DBH","HT","CANHT","ELEV","YEAR")

	dbh87 <- pp[trees.from.87,]$DBH87
	ht87 <- pp[trees.from.87,]$HT87
	elev87 <- pp[trees.from.87,]$ELEV
	plots87 <- pp[trees.from.87,]$PLOT
	canht87 <- pp[trees.from.87,]$CANHT87
	year87 <- rep(87,length(ht87))
	data87 <- data.frame(plots87,dbh87,ht87,canht87,elev87,year87)
	names(data87) <- c("PLOT","DBH","HT","CANHT","ELEV","YEAR")

	dbh98 <- pp[trees.from.98,]$DBH98
	ht98 <- pp[trees.from.98,]$HT98
	elev98 <- pp[trees.from.98,]$ELEV
	plots98 <- pp[trees.from.98,]$PLOT
	canht98 <- pp[trees.from.98,]$CANHT98
	year98 <- rep(98,length(ht98))
	data98 <- data.frame(plots98,dbh98,ht98,canht98,elev98,year98)
	names(data98) <- c("PLOT","DBH","HT","CANHT","ELEV","YEAR")

	dbh10 <- pp[trees.from.10,]$DBH10
	ht10 <- pp[trees.from.10,]$HT10
	elev10 <- pp[trees.from.10,]$ELEV
	plots10 <- pp[trees.from.10,]$PLOT
	canht10 <- pp[trees.from.10,]$CANHT10
	year10 <- rep(10,length(ht10))
	data10 <- data.frame(plots10,dbh10,ht10,canht10,elev10,year10)
	names(data10) <- c("PLOT","DBH","HT","CANHT","ELEV","YEAR")

	dataset <- rbind(data86,data87,data98,data10)
	dataset <- subset(dataset, PLOT > 3)
}

