calcpolinoms <- function(n) {
	res <- matrix(0.0, nrow = n + 1, ncol = n + 1)
	for(i in 0:n) {
		P <- Lpolinom(i)
		res[i+1,1:(i+1)] <- P
	}
	write.table(res, file = "Legendre_polinoms.txt", append = F, row.names = F, col.names = F)
	res
}




# function Lpolinom returns array of coefficients of Legendre's polinoms

Lpolinom <- function(deg) {
	P0 <- double()
	P1 <- double()
	P <- double()

	P0 <- 1.0			# degree = 0
	P1 <- c(0.0, 1.0)	# degree = 1

	if(deg == 0) {
		P0
	}
	else {

		if(deg == 1) {
			P1
		}
		else {			# for bigger degrees
				for(i in 2:deg) {	# calculations via recurence relation
					n <- i - 1
					P <- c(1:(i+1)) * 0.0
					P[2:(i+1)] <- P1[1:i] * (2 * n + 1) / (n + 1)
					P[1:(i-1)] <- P[1:(i-1)] - P0 * n / (n + 1)
					P0 <- P1
					P1 <- P
				}

			P1
		}
	}
}


# Pderivative returns array of coefficients of polinom
# which is m-th derivative of P

Pderivative <- function(P, m) {
	n <- length(P)		# n = degree of polinom + 1
	if(n <= m) {
		0
	}
	else {
		for(i in 1:m) {								# calculations of derivatives
			P <- P[2:(n-i+1)] * c(1:(n-i))	# step by step
		}
		P
	}
}


	# fact returns (l-m)! / (l+m)!

	fact <- function(l, m) {
		res <- 1.0
		if(m != 0) {
			for(i in (l+m):(l-m+1)) {
				res <- res / i
			}	
		}
		
		res
	}




# returns vallue of associated Legendre's function on x as argument
# Pl is a Legender's polinom degree = l

assLP <- function(Pl, l, m, x) {
	if(m == 0) {
		Pm <- Pl						# m-th derivative of Pl
	}
	else{ 
		Pm <- Pderivative(Pl, m)
	}
	Plmx <- sum(x^c(0:(l-m)) * Pm[1:(l-m+1)])		# Pm(x)
	Plm <- ((-1)^m) * Plmx * (1.0 - x * x)^(m/2.0)	# Plm(x)

	# we need normalization

	d0 <- 0.0 + (m == 0)
	Plm <- (-1.0)^m * sqrt((fact(l, m)) * (2.0 - d0)) * Plm

	if((abs(Plm) == Inf) | (sum(is.na(Plm) > 0))) {
		print(l, m, x)
		stop("Problem with Plm")
	}

	Plm
}

