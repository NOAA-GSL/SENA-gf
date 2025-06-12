AODC0 = 0.14  # Default value for aerosol optical depth
AODRETURN = 30.0  # Default value for AOD return time (minutes)
CCNCLEAN = max(5.0, (AODC0 / 0.0027) ** (1 / 0.640))
TCRIT = 258.0
TF = 258.16
TCR = 273.16  # Critical temperature for cloud water/ice conversion
TCRF = 1.0 / (TCR - TF)  # Scaling factor for temperature conversion
