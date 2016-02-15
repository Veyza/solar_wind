# Calculation of velosity of solar wind via the harmonic coefficients as input data

SWvelosity <- function(coefficients) 
{	
	R0 <- 2.5
	Rs <- 1.0
	r <- 0.5
	str_skip <- 17
	data <- read_harm_coef(coefficients, str_skip)
	g <- data$mat_g
	h <- data$mat_h

	P <- read_Lpolin(15)

	deg <- pi / 180.0			# 1 degree in radians
	theta <- c(1:180) * deg 	# latitude grid
	phi <- c(1:360) * deg 		# longitude grid

	  # Br <- read.csv("Radial_Field.csv")
	  # Br <- as.matrix(Br)
	  # Br <- Br[ ,-1]
	  Br <- outer(phi, theta, rfieldcalc, r, h, g, P)
	  write.table(Br, file = "Radial_Field.txt", append = F, col.names = F, row.names = F)
	  write.csv(Br,file = "Radial_Field.csv", col.names = F, row.names = F)

	  # Br <- rfieldcalc(phi[45], theta[135], r, h, g, P)
	  # Br

	  # Bth <- outer(phi, theta, thfieldcalc, r, h, g, P)
	  # write.table(Bth, file = "Latitude_Field.txt", append = F, col.names = F, row.names = F)
	  # write.csv(Bth, file = "Latitude_Field.csv", col.names = F, row.names = F)

	  # Bph <- outer(phi, theta, phfieldcalc, r, h, g, P)
	  # write.table(Bph, file = "Latitude_Field.txt", append = F, col.names = F, row.names = F)
	  # write.csv(Bph, file = "Latitude_Field.csv", col.names = F, row.names = F)

	  # Bl <- outer(phi, theta, lfieldcalc, r, h, g, P)

 
	  theta <- theta  / pi * 180
	  phi <- phi / pi * 180

	 png("./mapBr.png")
	 	par(fin = c(8, 4))
	 	filled.contour(phi, theta, Br, xlim = range(phi), ylim = range(theta))
	 dev.off()

	 #  png("./mapBth.png")
	 # 	par(fin = c(8, 4))
	 # 	filled.contour(phi, theta, Bth, xlim = range(phi), ylim = range(theta))
	 # dev.off()

	 # png("./mapBph.png")
	 # 	par(fin = c(8, 4))
	 # 	filled.contour(phi, theta, Bph, xlim = range(phi), ylim = range(theta))
	 # dev.off()

	 #  png("./mapBl.png")
	 # 	par(fin = c(8, 4))
	 # 	filled.contour(phi, theta, Bl, xlim = range(phi), ylim = range(theta))
	 # dev.off()

}


read_harm_coef <- function(input, str_skip) {
	g <- matrix(ncol = 16, nrow = 16)
	h <- matrix(ncol = 16, nrow = 16)
	for(i in 1 : 16) {
		data_g <- scan(file = input, what = double(), nmax = (i+1), skip = (i-1))
		g[i, 1:i] <- data_g[2:(i+1)]
		data_h <- scan(file = input, what = double(), nmax = (i+1), skip = (str_skip + i-1))
		h[i, 1:i] <- data_h[2:(i+1)]
	}
	output <- list(mat_g = g, mat_h = h)
	output
}


# n is the maximum degree of polinom
# n+1 polinoms will be read

read_Lpolin <- function(n) {
	input <- "Legendre_polinoms.txt"
	Lpolin <- matrix(0, nrow = n+1, ncol = n+1)
	for(i in 1:(n+1)) {
		x <- scan(file = input, what = double(), nmax = (n+1), skip = (i-1))
		Lpolin[i, ] <- x
	}
	Lpolin
}

rfieldcalc <- function(phi, theta, r, h, g, P) {
	R0 <- 2.5
	Rs <- 1.0
	
	Br <- 0.0
	mu <- cos(theta)
	Plm <- double()

	for(l in 0:15){
		item <- 0.0
		for(m in 0:l){
			Plm <- assLP(P[l+1,1:(l+1)], l, m, mu)
			item <- item + Plm * (g[l+1,m+1] * cos(phi * m) + h[l+1,m+1] * sin(phi * m))
				if(is.na(item)) {
					print(l)
					print(m)
					stop("current item is NA")
				}
		}
		Br <- Br + item * (R0 / r)^(l+2) * (l + 1 + l * (r / Rs)^(2*l+1)) / (l + 1 + l * (R0 / Rs)^(2*l+1))
	}

		# if(sum(abs(Br) == Inf) > 0) {
		# 	t <- theta[which(abs(Br) == Inf)]
		# 	ph <- phi[which(abs(Br) == Inf)]
		# }

	Br
}



thfieldcalc <- function(phi, theta, r, h, g, P) {
	R0 <- 2.5
	Rs <- 1.0
	
	Bth <- 0.0
	dtheta <- 1e-7
	mu <- cos(theta + dtheta)
	mu1 <- cos(theta - dtheta)
	dPlm <- double()

	for(l in 0:15){
	item <- 0.0
		for(m in 0:l){
			dPlm <- (assLP(P[l+1,1:(l+1)], l, m, mu) - assLP(P[l+1,1:(l+1)], l, m, mu1)) / (2.0 * dtheta)
			item <- item - dPlm * (g[l+1,m+1] * cos(phi * m) +h[l+1,m+1] * sin(phi * m))
		}
	Bth <- Bth + item * (R0 / r)^(l+2) * (1 - (r / Rs)^(2*l+1)) / (l + 1 + l * (R0 / Rs)^(2*l+1))
	}

	Bth

}



phfieldcalc <- function(phi, theta, r, h, g, P) {
	R0 <- 2.5
	Rs <- 1.0
	r <- R0

	Bph <- 0.0
	item <- 0.0
	mu <- cos(theta)

	for(l in 0:15) {
		item <- 0.0
		for(m in 0:l) {
			Plm <- assLP(P[l+1,1:(l+1)], l, m, mu)
			item <- item + m * Plm * (g[l+1,m+1] * sin(m * phi) - h[l+1,m+1] * cos(m * phi))
		}
		Bph <- Bph + item * (R0 / r)^(l+2) * (1.0 - (r / Rs)^(2*l+1)) / (l + 1 + l * (R0 / Rs)^(2*l+1))
	}

	Bph <- Bph / sin(theta)
	Bph
}


lfieldcalc <- function(phi, theta, r, h, g, P) {
	R0 <- 2.5
	Rs <- 1.0
	
	Br <- 0.0
	mu <- cos(theta)
	Plm <- double()

	for(l in 0:15){
		item <- 0.0
		for(m in 0:l){
			Plm <- assLP(P[l+1,1:(l+1)], l, m, mu)
			item <- item + Plm * (g[l+1,m+1] * cos(phi * m) +h[l+1,m+1] * sin(phi * m))
		}
		Br <- Br + item * (R0 / r)^(l+2) * (l + 1 + l * (r / Rs)^(2*l+1)) / (l + 1 + l * (R0 / Rs)^(2*l+1))
	}


	Bl <- Br * sin(theta) * cos(phi)

	Bl
}