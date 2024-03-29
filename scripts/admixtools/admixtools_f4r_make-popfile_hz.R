## Define output file
outfile_lines <- 'results/admixtools/input/popfile_f4r.txt'

## Function to create one input line
f4r <- function(A, B, X, C, O) {
  (pt1_3 <- paste(A, O))
  (pt2 <- paste(X, C))
  (pt4 <- paste(B, C))
  paste(pt1_3, ':', pt2, '::', pt1_3, ':', pt4)
}

## Create input lines one at a time
lines <- f4r(A = 'mmur_w', B = 'mmur_gan', X = 'mmur_hz', C = 'mgri_se', O = 'mruf')
lines <- c(lines, f4r(A = 'mmur_w', B = 'mmur_se', X = 'mmur_hz', C = 'mgri_se', O = 'mruf'))
lines <- c(lines, f4r(A = 'mmur_w', B = 'mmur_gan', X = 'mmur_hz', C = 'mgri_hz', O = 'mruf'))
lines <- c(lines, f4r(A = 'mmur_w', B = 'mmur_se', X = 'mmur_hz', C = 'mgri_hz', O = 'mruf'))
lines <- c(lines, f4r(A = 'mmur_w', B = 'mmur_gan', X = 'mmur_hz', C = 'mgri_sw', O = 'mruf'))
lines <- c(lines, f4r(A = 'mmur_w', B = 'mmur_se', X = 'mmur_hz', C = 'mgri_sw', O = 'mruf'))
lines <- c(lines, f4r(A = 'mmur_w', B = 'mmur_gan', X = 'mmur_se', C = 'mgri_se', O = 'mruf'))
lines <- c(lines, f4r(A = 'mmur_w', B = 'mmur_gan', X = 'mmur_se', C = 'mgri_hz', O = 'mruf'))
lines <- c(lines, f4r(A = 'mmur_w', B = 'mmur_gan', X = 'mmur_se', C = 'mgri_sw', O = 'mruf'))
lines <- c(lines, f4r(A = 'mgri_sw', B = 'mgri_hz', X = 'mmur_gan', C = 'mmur_w', O = 'mruf'))
lines <- c(lines, f4r(A = 'mgri_sw', B = 'mgri_hz', X = 'mmur_se', C = 'mmur_w', O = 'mruf'))
lines <- c(lines, f4r(A = 'mgri_sw', B = 'mgri_hz', X = 'mmur_se', C = 'mmur_gan', O = 'mruf'))
lines <- c(lines, f4r(A = 'mgri_sw', B = 'mgri_se', X = 'mmur_gan', C = 'mmur_w', O = 'mruf'))
lines <- c(lines, f4r(A = 'mgri_sw', B = 'mgri_se', X = 'mmur_se', C = 'mmur_w', O = 'mruf'))
lines <- c(lines, f4r(A = 'mgri_sw', B = 'mgri_se', X = 'mmur_se', C = 'mmur_gan', O = 'mruf'))

## Write to file
writeLines(lines, outfile_lines)
