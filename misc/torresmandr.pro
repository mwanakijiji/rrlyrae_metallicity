PRO torresmandr, t_eff, t_eff_err, logg, logg_err, feh, feherr
  ;;applies Torres et al. 2010 polynomial fit given T_eff, log(g) and [Fe/H] to obtain mass and radius of primary, based on observations of eclipsing binaries
  
  x = DOUBLE(ALOG10(t_eff)-4.1D) ;;as in Torres paper
  xerr = DOUBLE(1.D/ALOG(10.D)*(t_eff_err/t_eff))

  ;;coefficients from Torres paper
  a1 = 1.5689D & a2 = 1.3787D & a3 = 0.4243D & a4 = 1.139D & a5 = -0.1425D & a6 = 0.01969D & a7 = 0.1010D
  b1 = 2.4427D & b2 = 0.6679D & b3 = 0.1771D & b4 = 0.705D & b5 = -0.21415D & b6 = 0.02306D & b7 = 0.04173D

  torreslogmscatter = 0.027 & torreslogrscatter = 0.014 ;;scatter reported in Torres et al. paper, to be added in quadrature to rest of error at end

  ;;read in covariance matrix for logM
  OPENR, flun, "torres_logM_covar.txt", /GET_LUN & logMcovar = DBLARR(7,7)
  READF, flun, logMcovar
  FREE_LUN, flun
  ;;Willie Torres scaled the original data by 1000 to preserve sig. digits, so convert back
  logMcovar = logMcovar/(1000.D)

  ;;read in covariance matrix for logR
  OPENR, flun, "torres_logR_covar.txt", /GET_LUN & logRcovar = DBLARR(7,7)
  READF, flun, logRcovar
  FREE_LUN, flun
  ;;Willie Torres scaled the original data by 1000 to preserve sig. digits, so convert back
  logRcovar = logRcovar/(1000.D)

  ;;calculate log(M)
  logm = a1 + a2*x + a3*x^2 + a4*x^3 + a5*(logg)^2 + a6*(logg)^3 + a7*feh
  ;;calculate log(R)
  logr = b1 + b2*x + b3*x^2 + b4*x^3 + b5*(logg)^2 + b6*(logg)^3 + b7*feh

  ;;calculate covariance terms for adding in final equation for error in log(M)
  covarerr_a12 = 2.D * (1.D) * (x) * (logMcovar[0,1])
  covarerr_a13 = 2.D * (1.D) * (x^2) * (logMcovar[0,2])
  covarerr_a14 = 2.D * (1.D) * (x^3) * (logMcovar[0,3])
  covarerr_a15 = 2.D * (1.D) * (logg^2) * (logMcovar[0,4])
  covarerr_a16 = 2.D * (1.D) * (logg^3) * (logMcovar[0,5])
  covarerr_a17 = 2.D * (1.D) * (feh) * (logMcovar[0,6])
  covarerr_a23 = 2.D * (x) * (x^2) * (logMcovar[1,2])
  covarerr_a24 = 2.D * (x) * (x^3) * (logMcovar[1,3])
  covarerr_a25 = 2.D * (x) * (logg^2) * (logMcovar[1,4])
  covarerr_a26 = 2.D * (x) * (logg^3) * (logMcovar[1,5])
  covarerr_a27 = 2.D * (x) * (feh) * (logMcovar[1,6])
  covarerr_a34 = 2.D * (x^2) * (x^3) * (logMcovar[2,3])
  covarerr_a35 = 2.D * (x^2) * (logg^2) * (logMcovar[2,4])
  covarerr_a36 = 2.D * (x^2) * (logg^3) * (logMcovar[2,5])
  covarerr_a37 = 2.D * (x^2) * (feh) * (logMcovar[2,6])
  covarerr_a45 = 2.D * (x^3) * (logg^2) * (logMcovar[3,4])
  covarerr_a46 = 2.D * (x^3) * (logg^3) * (logMcovar[3,5])
  covarerr_a47 = 2.D * (x^3) * (feh) * (logMcovar[3,6])
  covarerr_a56 = 2.D * (logg^2) * (logg^3) * (logMcovar[4,5])
  covarerr_a57 = 2.D * (logg^2) * (feh) * (logMcovar[4,6])
  covarerr_a67 = 2.D * (logg^3) * (feh) * (logMcovar[5,6])

  ;;calculate error in log(M)
  errlogm = SQRT( (logMcovar[0,0])*(1.)^2 + (logMcovar[1,1])*(x)^2 + (logMcovar[2,2])*(x^2)^2 + (logMcovar[3,3])*(x^3)^2 + (logMcovar[4,4])*(logg^2)^2 + (logMcovar[5,5])*(logg^3)^2 + (logMcovar[6,6])*(feh)^2 + (xerr)^2*(a2+2.*a3*x+3*a4*x^2)^2 + (logg_err)^2*(2.*a5*logg+3.*a6*logg^2)^2 + (feherr)^2*(a7)^2 + covarerr_a12 + covarerr_a13 + covarerr_a14 + covarerr_a15 + covarerr_a16 + covarerr_a17 + covarerr_a23 + covarerr_a24 + covarerr_a25 + covarerr_a26 + covarerr_a27 + covarerr_a34 + covarerr_a35 + covarerr_a36 + covarerr_a37 + covarerr_a45 + covarerr_a46 + covarerr_a47 + covarerr_a56 + covarerr_a57 + covarerr_a67)

  ;;add scatter from Torres paper in quadrature
  errlogm = SQRT(errlogm^2. + torreslogmscatter^2.)

  ;;calculate covariance terms for adding in final equation for error in log(R)
  covarerr_b12 = 2.D * (1.D) * (x) * (logRcovar[0,1])
  covarerr_b13 = 2.D * (1.D) * (x^2) * (logRcovar[0,2])
  covarerr_b14 = 2.D * (1.D) * (x^3) * (logRcovar[0,3])
  covarerr_b15 = 2.D * (1.D) * (logg^2) * (logRcovar[0,4])
  covarerr_b16 = 2.D * (1.D) * (logg^3) * (logRcovar[0,5])
  covarerr_b17 = 2.D * (1.D) * (feh) * (logRcovar[0,6])
  covarerr_b23 = 2.D * (x) * (x^2) * (logRcovar[1,2])
  covarerr_b24 = 2.D * (x) * (x^3) * (logRcovar[1,3])
  covarerr_b25 = 2.D * (x) * (logg^2) * (logRcovar[1,4])
  covarerr_b26 = 2.D * (x) * (logg^3) * (logRcovar[1,5])
  covarerr_b27 = 2.D * (x) * (feh) * (logRcovar[1,6])
  covarerr_b34 = 2.D * (x^2) * (x^3) * (logRcovar[2,3])
  covarerr_b35 = 2.D * (x^2) * (logg^2) * (logRcovar[2,4])
  covarerr_b36 = 2.D * (x^2) * (logg^3) * (logRcovar[2,5])
  covarerr_b37 = 2.D * (x^2) * (feh) * (logRcovar[2,6])
  covarerr_b45 = 2.D * (x^3) * (logg^2) * (logRcovar[3,4])
  covarerr_b46 = 2.D * (x^3) * (logg^3) * (logRcovar[3,5])
  covarerr_b47 = 2.D * (x^3) * (feh) * (logRcovar[3,6])
  covarerr_b56 = 2.D * (logg^2) * (logg^3) * (logRcovar[4,5])
  covarerr_b57 = 2.D * (logg^2) * (feh) * (logRcovar[4,6])
  covarerr_b67 = 2.D * (logg^3) * (feh) * (logRcovar[5,6])

  ;;calculate error in log(R)
  errlogr = SQRT( (logRcovar[0,0])*(1.)^2 + (logRcovar[1,1])*(x)^2 + (logRcovar[2,2])*(x^2)^2 + (logRcovar[3,3])*(x^3)^2 + (logRcovar[4,4])*(logg^2)^2 + (logRcovar[5,5])*(logg^3)^2 + (logRcovar[6,6])*(feh)^2 + (xerr)^2*(b2+2.*b3*x+3*b4*x^2)^2 + (logg_err)^2*(2.*b5*logg+3.*b6*logg^2)^2 + (feherr)^2*(b7)^2 + covarerr_b12 + covarerr_b13 + covarerr_b14 + covarerr_b15 + covarerr_b16 + covarerr_b17 + covarerr_b23 + covarerr_b24 + covarerr_b25 + covarerr_b26 + covarerr_b27 + covarerr_b34 + covarerr_b35 + covarerr_b36 + covarerr_b37 + covarerr_b45 + covarerr_b46 + covarerr_b47 + covarerr_b56 + covarerr_b57 + covarerr_b67 )

  ;;add scatter from Torres paper in quadrature
  errlogr = SQRT(errlogr^2. + torreslogrscatter^2.)

  print, "Independent error terms for logM are = "+STRTRIM(STRING((logMcovar[0,0])*(1.)^2 + (logMcovar[1,1])*(x)^2 + (logMcovar[2,2])*(x^2)^2 + (logMcovar[3,3])*(x^3)^2 + (logMcovar[4,4])*(logg^2)^2 + (logMcovar[5,5])*(logg^3)^2 + (logMcovar[6,6])*(feh)^2 + (xerr)^2*(a2+2.*a3*x+3*a4*x^2)^2 + (logg_err)^2*(2.*a5*logg+3.*a6*logg^2)^2 + (feherr)^2*(a7)^2,F='(D)'),2)
  print, "Dependent error terms for logM are = "+STRTRIM(STRING(covarerr_a12 + covarerr_a13 + covarerr_a14 + covarerr_a15 + covarerr_a16 + covarerr_a17 + covarerr_a23 + covarerr_a24 + covarerr_a25 + covarerr_a26 + covarerr_a27 + covarerr_a34 + covarerr_a35 + covarerr_a36 + covarerr_a37 + covarerr_a45 + covarerr_a46 + covarerr_a47 + covarerr_a56 + covarerr_a57 + covarerr_a67,F='(D)'),2)
  print, "Total error in logM is = "+STRTRIM(STRING(errlogm,F='(D)'),2)
  
  print, "Independent error terms for logR are = "+STRTRIM(STRING((logRcovar[0,0])*(1.)^2 + (logRcovar[1,1])*(x)^2 + (logRcovar[2,2])*(x^2)^2 + (logRcovar[3,3])*(x^3)^2 + (logRcovar[4,4])*(logg^2)^2 + (logRcovar[5,5])*(logg^3)^2 + (logRcovar[6,6])*(feh)^2 + (xerr)^2*(b2+2.*b3*x+3*b4*x^2)^2 + (logg_err)^2*(2.*b5*logg+3.*b6*logg^2)^2 + (feherr)^2*(b7)^2,F='(D)'),2)
  print, "Dependent error terms for logR are = "+STRTRIM(STRING(covarerr_b12 + covarerr_b13 + covarerr_b14 + covarerr_b15 + covarerr_b16 + covarerr_b17 + covarerr_b23 + covarerr_b24 + covarerr_b25 + covarerr_b26 + covarerr_b27 + covarerr_b34 + covarerr_b35 + covarerr_b36 + covarerr_b37 + covarerr_b45 + covarerr_b46 + covarerr_b47 + covarerr_b56 + covarerr_b57 + covarerr_b67,F='(D)'),2)
  print, "Total error in logR is = "+STRTRIM(STRING(errlogr,F='(D)'),2)

  print, "log(mass) is "+STRTRIM(STRING(logm,F='(F9.6)'),2)+" +/- "+STRTRIM(STRING(errlogm,F='(F8.6)'),2)+"."
  print, "log(radius) is "+STRTRIM(STRING(logr,F='(F9.6)'),2)+" +/- "+STRTRIM(STRING(errlogr,F='(F8.6)'),2)+"."

  print, "The mass is "+STRTRIM(STRING(10.^(logm)),2)+" +/- "+STRTRIM(STRING((10.^(logm+errlogm) - 10^(logm-errlogm))/2.),2)+" solar masses."
  print, "The radius is "+STRTRIM(STRING(10.^(logr)),2)+" +/- "+STRTRIM(STRING((10.^(logr+errlogr) - 10.^(logr-errlogr))/2.),2)+" solar radii."

;  print, "The mass is "+STRTRIM(STRING(10.^(logm-errlogm)),2)+", "+STRTRIM(STRING(10.^(logm)),2)+", "+STRTRIM(STRING(10.^(logm+errlogm)),2)+" solar masses."
;  print, "The radius is "+STRTRIM(STRING(10.^(logr-errlogr)),2)+", "+STRTRIM(STRING(10.^(logr)),2)+", "+STRTRIM(STRING(10.^(logr+errlogr)),2)+" solar radii."
  RETURN
END
