











# generate the harmonic coefficients representing the field of the dipole
# q is an absolute value of point charges, dl is a distanse between the charges
# r is a distanse from dipole on which we need to compute field

coefDfield <- function(q = 1, dl = 10, r = 1e+3) {

	d <- q * dl
	mu <- 4.0 * pi * 1e-7		# magnetic permeability of vacuum (SI)
	theta <- c(1:180) * pi / 180.0
	phi <- c(1:360) * pi / 180.0
		B <- as.matrix(read.csv("Radial_Field.csv"))


		harmCoef <- function(B, l, m, theta, phi, r, d) {
			Pl <- Lpolinom(l)
			n <- length(B[1, ])
			k <- length(B[ ,1])
			print(n)
			print(k)
			if((n == 360 | k == 360) && (n == 180 | k == 180)) {
				X <- 180
				Y <- 360
				gml <- 0.0
				hml <- 0.0
					for(i in 1:X) {
						for(j in 1:Y) {
							# B <-d / mu / r^3 * sqrt(1.0 + 3.0 * cos(theta[i])^2)
							if(is.na(B[j,i])) {
								print(i)
								print(j)
							}
							gml <- gml + B[j,i] * assLP(Pl, l, m, cos(theta[i])) * (cos(m * phi[j]))
							hml <- hml + B[j,i] * assLP(Pl, l, m, cos(theta[i])) * (sin(m * phi[j]))
						}
					}
				gml <- gml * (2.0 * l) / X / Y
				hml <- hml * (2.0 * l) / X / Y
				c(gml, hml)
		}
		else{
			print("X______X")
		}

		}



	for(l in 0:15) {
		g <- double()
		h <- double()
		for(m in 0:l) {
			gh <- harmCoef(B, l, m, theta, phi, r, d)
			g[m+1] <- gh[1]
			h[m+1] <- gh[2]
		}
		print(l)
		print(m)
		 write(c(l, g), "dipole_g.grm", append = T, sep = "  ", ncolumns = 17)
		 write(c(l, h), "dipole_h.grm", append = T, sep = "  ", ncolumns = 17)
	}




}