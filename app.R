# =============================================================================
# Time-to-Event Sample Size Calculator — R Shiny translation of .jsx file
# Exponential · Log-rank · Schoenfeld / Lachin-Foulkes
# Author: Mary Ann Binuya
# Date v1: 2025-10-08
# --one-arm
# Date v2: 2026-03-26
# --update: added two-arm (formulas WIP/not final; some tests are inherently one-sided)
# =============================================================================
# Required packages:
#   shiny, shinyjs, htmltools
# Install with: install.packages(c("shiny","shinyjs"))
# =============================================================================

library(shiny)
library(shinyjs)

# ─── Colour palette (mirroring JSX constants) ─────────────────────────────────
BG  <- "#030712"; PAN <- "#0a1628"; BDR <- "#0f2137"; B2  <- "#1e3a5f"
TX  <- "#e2e8f0"; MU  <- "#64748b"; DM  <- "#334155"
BL  <- "#60a5fa"; PK  <- "#f472b6"; GN  <- "#34d399"; OR  <- "#fb923c"
GY  <- "#94a3b8"; RD  <- "#f87171"; TL  <- "#2dd4bf"

# ─── Math helpers ─────────────────────────────────────────────────────────────
# R's built-in qnorm / dnorm / pnorm are used directly.
# They are more accurate than the JS Wichura polynomial — functionally identical.

lm_fn  <- function(m)   if (m > 0)  log(2) / m  else NaN   # median -> lambda
ml_fn  <- function(l)   if (l > 0)  log(2) / l  else NaN   # lambda -> median
ls_fn  <- function(S,t) if (S>0 && S<1 && t>0) -log(S)/t else NaN  # S(t*) -> lambda
Sv_fn  <- function(l,t) exp(-l*t)                            # S(t)

# P(event) under uniform accrual A over [0,A], follow-up to L
pEv_fn <- function(lam, A, L) {
  if (!is.finite(lam) || lam <= 0) return(NaN)
  F_ <- L - A
  if (A <= 0) return(max(0, min(1, 1 - exp(-lam * L))))
  max(0, min(1,
             1 - (exp(-lam * F_) - exp(-lam * L)) / (A * lam)
  ))
}

# Cumulative event probability at calendar time t (Lachin-Foulkes)
cumP_fn <- function(t, lam, A) {
  if (!is.finite(lam) || lam <= 0 || t <= 0) return(0)
  if (A <= 0) return(min(1, 1 - exp(-lam * t)))
  if (t <= A) return(max(0, 1 - (1 - exp(-lam * t)) / (lam * t)))
  max(0,
      1 - (exp(-lam * (t - A)) - exp(-lam * t)) / (A * lam)
  )
}

# ─── Single-arm calculation ────────────────────────────────────────────────────
calcOne <- function(l0, l1, al, pw, sd, A, L, dr) {
  HR <- l1 / l0
  if (!is.finite(HR) || abs(log(HR)) < 1e-9) return(NULL)
  z1a <- if (sd == 2) qnorm(1 - al / 2) else qnorm(1 - al)
  z1b <- qnorm(pw)
  lHR <- log(HR)
  D   <- ceiling(((z1a + z1b) / lHR)^2)
  pe  <- pEv_fn(l1, A, L)
  peff<- pe * (1 - dr / 100)
  N   <- if (peff > 0) ceiling(D / peff) else NaN
  mu  <- sqrt(D) * abs(lHR)
  apw <- if (sd == 2) pnorm(mu - z1a) + pnorm(-mu - z1a) else pnorm(mu - z1a)
  list(HR=HR, D=D, pe=pe, peff=peff, N=N, z1a=z1a, z1b=z1b,
       lHR=lHR, mu=mu, pw=apw,
       pc=NULL, pt=NULL, pAv=NULL, Nt=NULL, nt=NULL, nc=NULL, eLHR=NULL)
}

# ─── Two-arm calculation ───────────────────────────────────────────────────────
calcTwo <- function(l0, l1, al, pw, sd, obj, dlt, k, A, L, dr) {
  HR  <- l1 / l0
  if (!is.finite(HR) || HR <= 0) return(NULL)
  lHR <- log(HR)
  r   <- k / (1 + k)
  
  if (obj == "equivalence") {
    ef  <- abs(dlt) - abs(lHR)
    if (ef <= 0) return(NULL)
    za  <- qnorm(1 - al)
    zb  <- qnorm(pw)
    D   <- ceiling((za + zb)^2 / (ef^2 * r * (1 - r)))
    pc  <- pEv_fn(l0, A, L); pt <- pEv_fn(l1, A, L)
    pAv <- r * pt + (1 - r) * pc
    peff<- pAv * (1 - dr / 100)
    Nt  <- if (peff > 0) ceiling(D / peff) else NaN
    nt  <- ceiling(Nt * r); nc <- Nt - nt
    mu  <- sqrt(D) * ef * sqrt(r * (1 - r))
    apw <- pnorm(mu - za) - pnorm(-mu - za)
    return(list(HR=HR, D=D, pc=pc, pt=pt, pAv=pAv, peff=peff,
                Nt=Nt, nt=nt, nc=nc,
                z1a=za, z1b=zb, lHR=lHR, eLHR=ef, mu=mu, pw=apw,
                pe=NULL, N=NULL))
  }
  
  eL <- switch(obj,
               "equality"       = lHR,
               "noninferiority" = dlt - lHR,
               "superiority"    = -(lHR + dlt),
               return(NULL)
  )
  z1a <- if (obj == "equality" && sd == 2) qnorm(1 - al / 2) else qnorm(1 - al)
  if (abs(eL) < 1e-9) return(NULL)
  
  z1b <- qnorm(pw)
  D   <- ceiling((z1a + z1b)^2 / (eL^2 * r * (1 - r)))
  pc  <- pEv_fn(l0, A, L); pt <- pEv_fn(l1, A, L)
  pAv <- r * pt + (1 - r) * pc
  peff<- pAv * (1 - dr / 100)
  Nt  <- if (peff > 0) ceiling(D / peff) else NaN
  nt  <- ceiling(Nt * r); nc <- Nt - nt
  mu  <- sqrt(D) * abs(eL) * sqrt(r * (1 - r))
  apw <- if (sd == 2) pnorm(mu - z1a) + pnorm(-mu - z1a) else pnorm(mu - z1a)
  list(HR=HR, D=D, pc=pc, pt=pt, pAv=pAv, peff=peff,
       Nt=Nt, nt=nt, nc=nc,
       z1a=z1a, z1b=z1b, lHR=lHR, eLHR=eL, mu=mu, pw=apw,
       pe=NULL, N=NULL)
}

# ─── Budget back-calculation ───────────────────────────────────────────────────
backCalc <- function(Nbudget, l0, l1, al, pw, sd, obj, dlt, k, A, L, dr, arms) {
  if (!is.finite(Nbudget) || Nbudget < 1) return(NULL)
  HR <- l1 / l0
  if (!is.finite(HR) || HR <= 0 || abs(log(HR)) < 1e-9) return(NULL)
  lHR <- log(HR)
  
  effLHR <- if (arms == 1) {
    abs(lHR)
  } else {
    switch(obj,
           "equality"       = abs(lHR),
           "noninferiority" = abs(dlt - lHR),
           "superiority"    = abs(lHR + dlt),
           "equivalence"    = max(0, abs(dlt) - abs(lHR)),
           return(NULL)
    )
  }
  if (effLHR < 1e-9) return(NULL)
  
  pe <- if (arms == 1) {
    pEv_fn(l1, A, L) * (1 - dr / 100)
  } else {
    r <- k / (1 + k)
    (r * pEv_fn(l1, A, L) + (1 - r) * pEv_fn(l0, A, L)) * (1 - dr / 100)
  }
  if (!is.finite(pe) || pe <= 0) {
    return(list(feasible=FALSE, reason="p(event) is zero — extend follow-up or adjust hazards."))
  }
  
  Nmin <- ceiling(1 / pe)
  if (Nbudget < Nmin) {
    return(list(
      feasible=FALSE,
      reason=paste0("With p\u2091 = ", sprintf("%.1f%%", pe*100),
                    ", fewer than 1 event is expected at N = ", Nbudget,
                    ". Minimum feasible N is ", Nmin, ".")
    ))
  }
  
  D  <- floor(Nbudget * pe)
  r2 <- if (arms == 2) k / (1 + k) else 1
  mu1 <- if (arms == 1) sqrt(D) * effLHR else sqrt(D) * effLHR * sqrt(r2 * (1 - r2))
  
  z1a_orig <- if (sd == 2) qnorm(1 - al / 2) else qnorm(1 - al)
  power    <- if (sd == 2) pnorm(mu1 - z1a_orig) + pnorm(-mu1 - z1a_orig)
  else         pnorm(mu1 - z1a_orig)
  
  z1b        <- qnorm(pw)
  z1a_needed <- mu1 - z1b
  alpha1s    <- if (z1a_needed > 0) 1 - pnorm(z1a_needed) else NULL
  alpha2s    <- if (z1a_needed > 0) 2 * (1 - pnorm(z1a_needed)) else NULL
  impossible <- z1a_needed <= 0
  
  list(feasible=TRUE, D=D, mu1=mu1, power=power,
       alpha1s=alpha1s, alpha2s=alpha2s,
       z1a_needed=z1a_needed, impossible=impossible, pe=pe)
}

# ─── Power given N (for heatmap) ──────────────────────────────────────────────
pwrN_fn <- function(N, HR, al, pe, sd, r=NULL) {
  if (!is.finite(pe) || pe <= 0 || !is.finite(HR) || HR <= 0) return(0)
  if (abs(log(HR)) < 1e-9) return(0)
  z      <- if (sd == 2) qnorm(1 - al / 2) else qnorm(1 - al)
  D      <- N * pe
  rFactor<- if (!is.null(r) && r > 0 && r < 1) sqrt(r * (1 - r)) else 1
  mu     <- sqrt(D) * abs(log(HR)) * rFactor
  if (sd == 2) pnorm(mu - z) + pnorm(-mu - z) else pnorm(mu - z)
}

# ─── SVG builders ─────────────────────────────────────────────────────────────

make_surv_svg <- function(l0, l1, A, L, mode, ts, arms, unit_short, unit_long) {
  W=520; H=280; PL=50; PR=14; PT=36; PB=42
  w <- W - PL - PR; h <- H - PT - PB
  m0 <- ml_fn(l0); m1 <- ml_fn(l1)
  tMax <- max(L + 4,
              if (is.finite(m0)) m0 * 3.2 else 0,
              if (is.finite(m1)) m1 * 3.2 else 0, 24)
  xS <- function(t) PL + (t / tMax) * w
  yS <- function(s) PT + h - s * h
  mkPath <- function(lam) {
    pts <- sapply(seq(0, 160) / 160 * tMax, function(t)
      paste0(sprintf("%.1f", xS(t)), ",", sprintf("%.1f", yS(Sv_fn(lam, t)))))
    paste0("M", paste(pts, collapse="L"))
  }
  xTk <- round(seq(0, 6) / 6 * tMax)
  yTk <- c(0, 0.25, 0.5, 0.75, 1.0)
  
  s <- sprintf('<svg viewBox="0 0 %d %d" style="width:100%%;height:auto">', W, H)
  
  for (y in yTk)
    s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>',
                           PL, PL+w, yS(y), yS(y), BDR))
  for (x in xTk)
    s <- paste0(s, sprintf('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="1"/>',
                           xS(x), xS(x), PT, PT+h, BDR))
  if (A > 0) {
    aw <- xS(min(A, tMax)) - xS(0)
    s <- paste0(s, sprintf('<rect x="%.1f" y="%d" width="%.1f" height="%d" fill="%s10"/>',
                           xS(0), PT, aw, h, PK))
  }
  if (L <= tMax)
    s <- paste0(s, sprintf('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="1" stroke-dasharray="3,3"/>',
                           xS(L), xS(L), PT, PT+h, DM))
  s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="%s" stroke-width="1.5"/>',
                         PL, PL+w, PT+h, PT+h, GY))
  s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="%s" stroke-width="1.5"/>',
                         PL, PL, PT, PT+h, GY))
  s <- paste0(s, sprintf('<path d="%s" stroke="%s" stroke-width="2" fill="none"/>', mkPath(l0), GY))
  s <- paste0(s, sprintf('<path d="%s" stroke="%s" stroke-width="2" fill="none"/>', mkPath(l1), OR))
  # median line at S=0.5
  s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="%s" stroke-width="1" stroke-dasharray="5,3"/>',
                         PL, PL+w, yS(0.5), yS(0.5), DM))
  # median drop-lines
  for (pair in list(list(lam=l0, col=GY, up=TRUE), list(lam=l1, col=OR, up=FALSE))) {
    m <- ml_fn(pair$lam)
    if (!is.finite(m) || m > tMax) next
    s <- paste0(s, sprintf('<line x1="%.1f" x2="%.1f" y1="%.1f" y2="%.1f" stroke="%s" stroke-width="1" stroke-dasharray="4,3"/>',
                           xS(m), xS(m), yS(0.5), yS(0), pair$col))
    ry <- if (pair$up) yS(0.5) - 23 else yS(0.5) + 5
    ty <- if (pair$up) yS(0.5) - 11 else yS(0.5) + 17
    s <- paste0(s, sprintf('<rect x="%.1f" y="%.1f" width="64" height="16" rx="4" fill="%s" stroke="%s"/>',
                           xS(m)-32, ry, PAN, if(pair$up) DM else B2))
    s <- paste0(s, sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="9" fill="%s">m=%.1f %s</text>',
                           xS(m), ty, pair$col, m, unit_short))
  }
  if (mode == "St" && is.finite(ts) && ts > 0 && ts <= tMax) {
    s <- paste0(s, sprintf('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="1.2" stroke-dasharray="4,3"/>',
                           xS(ts), xS(ts), PT, PT+h, BL))
    s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="9" fill="%s">t*=%.1f %s</text>',
                           xS(ts), PT+15, BL, ts, unit_short))
    s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="8" fill="%s">S0=%.3f</text>',
                           xS(ts)-10, PT+27, GY, Sv_fn(l0, ts)))
    s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="8" fill="%s">S1=%.3f</text>',
                           xS(ts)+30, PT+27, OR, Sv_fn(l1, ts)))
  }
  for (x in xTk)
    s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="9" fill="%s">%g</text>',
                           xS(x), PT+h+14, MU, x))
  for (y in yTk)
    s <- paste0(s, sprintf('<text x="%d" y="%.1f" text-anchor="end" font-size="9" fill="%s">%g%%</text>',
                           PL-5, yS(y)+4, MU, y*100))
  s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="9" fill="%s">Time (%s)</text>',
                         PL+w/2, H-3, MU, unit_long))
  s <- paste0(s, sprintf('<text x="12" y="%.1f" text-anchor="middle" font-size="9" fill="%s" transform="rotate(-90,12,%.1f)">S(t)</text>',
                         PT+h/2, MU, PT+h/2))
  lb0 <- if (arms == 1) "H\u2080 Benchmark" else "Control"
  lb1 <- if (arms == 1) "H\u2081 Target"    else "Treatment"
  s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="%s" stroke-width="2"/>',
                         PL, PL+14, PT-14, PT-14, GY))
  s <- paste0(s, sprintf('<text x="%d" y="%d" font-size="9" fill="%s">%s</text>',
                         PL+17, PT-10, GY, lb0))
  s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="%s" stroke-width="2"/>',
                         PL+110, PL+124, PT-14, PT-14, OR))
  s <- paste0(s, sprintf('<text x="%d" y="%d" font-size="9" fill="%s">%s</text>',
                         PL+127, PT-10, OR, lb1))
  paste0(s, "</svg>")
}

make_z_svg <- function(mu1, z1a, sd) {
  if (!is.finite(mu1) || mu1 <= 0) return("")
  W=520; H=260; PL=50; PR=14; PT=36; PB=42
  w <- W-PL-PR; h <- H-PT-PB
  xMin <- min(-4, -z1a - 1)
  xMax <- max(z1a + 1, mu1 + 3.5)
  xR   <- xMax - xMin
  xS   <- function(x) PL + (x - xMin) / xR * w
  N    <- 260
  xs   <- seq(xMin, xMax, length.out=N)
  y0   <- dnorm(xs)
  y1   <- dnorm(xs, mu1)
  yMax <- max(y0, y1) * 1.2
  yS   <- function(y) PT + h - (y / yMax) * h
  base <- PT + h
  pw   <- if (sd == 2) pnorm(mu1-z1a)+pnorm(-mu1-z1a) else pnorm(mu1-z1a)
  
  mkFill <- function(ys, fn) {
    pxs  <- xS(xs)
    keep <- fn(pxs)
    seg  <- which(keep)
    if (length(seg) < 2) return("")
    pts <- paste0(sprintf("%.1f", pxs[seg]), ",", sprintf("%.1f", yS(ys[seg])))
    paste0("M", sprintf("%.1f", pxs[seg[1]]), ",", base,
           "L", paste(pts, collapse="L"),
           "L", sprintf("%.1f", pxs[seg[length(seg)]]), ",", base, "Z")
  }
  mkLine <- function(ys)
    paste0("M", paste0(sprintf("%.1f", xS(xs)), ",", sprintf("%.1f", yS(ys)), collapse="L"))
  
  xc <- xS(z1a); xcN <- xS(-z1a)
  s <- sprintf('<svg viewBox="0 0 %d %d" style="width:100%%;height:auto">', W, H)
  s <- paste0(s, sprintf('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="1" stroke-dasharray="2,3"/>',
                         xS(0), xS(0), PT, PT+h, DM))
  s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="%s" stroke-width="1.5"/>',
                         PL, PL+w, PT+h, PT+h, GY))
  s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="%s" stroke-width="1.5"/>',
                         PL, PL, PT, PT+h, GY))
  
  aR <- mkFill(y0, function(px) px >= xc)
  if (nchar(aR) > 0) s <- paste0(s, sprintf('<path d="%s" fill="%s77"/>', aR, GY))
  if (sd == 2) {
    aL <- mkFill(y0, function(px) px <= xcN)
    if (nchar(aL) > 0) s <- paste0(s, sprintf('<path d="%s" fill="%s77"/>', aL, GY))
  }
  bP <- if (sd == 2)
    mkFill(y1, function(px) px >= xcN & px <= xc)
  else
    mkFill(y1, function(px) px <= xc)
  if (nchar(bP) > 0) s <- paste0(s, sprintf('<path d="%s" fill="%s55"/>', bP, BL))
  
  s <- paste0(s, sprintf('<path d="%s" stroke="%s" stroke-width="2" fill="none"/>', mkLine(y0), GY))
  s <- paste0(s, sprintf('<path d="%s" stroke="%s" stroke-width="2" fill="none"/>', mkLine(y1), BL))
  s <- paste0(s, sprintf('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="1.5" stroke-dasharray="5,3"/>',
                         xc, xc, PT, PT+h, RD))
  s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="9" fill="%s">z=%.3f</text>',
                         xc, PT+12, RD, z1a))
  if (sd == 2) {
    s <- paste0(s, sprintf('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="1.5" stroke-dasharray="5,3"/>',
                           xcN, xcN, PT, PT+h, RD))
    s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="9" fill="%s">-z</text>',
                           xcN, PT+12, RD))
  }
  if (xS(mu1) < PL + w) {
    s <- paste0(s, sprintf('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="2"/>',
                           xS(mu1), xS(mu1), PT+h-4, PT+h+4, BL))
    s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="9" fill="%s">mu=%.2f</text>',
                           xS(mu1), PT+h+17, BL, mu1))
  }
  xTk <- seq(xMin, xMax, length.out=7)
  for (x in xTk)
    s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="8" fill="%s">%.1f</text>',
                           xS(x), PT+h+14, MU, x))
  s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="9" fill="%s">Z statistic</text>',
                         PL+w/2, H-3, MU))
  s <- paste0(s, sprintf('<rect x="%d" y="%d" width="78" height="38" rx="5" fill="%s" stroke="%s"/>',
                         W-PR-80, PT+4, PAN, B2))
  s <- paste0(s, sprintf('<text x="%d" y="%d" text-anchor="middle" font-size="9" fill="%s">Power=%.1f%%</text>',
                         W-PR-41, PT+18, GN, pw*100))
  s <- paste0(s, sprintf('<text x="%d" y="%d" text-anchor="middle" font-size="9" fill="%s">beta=%.1f%%</text>',
                         W-PR-41, PT+31, BL, (1-pw)*100))
  paste0(s, "</svg>")
}

make_evt_svg <- function(l1, A, L, dr, N, l0, arms, k, unit_short, unit_long) {
  if (!is.finite(l1) || !is.finite(N) || N <= 0 || L <= 0) return("")
  W=520; H=240; PL=58; PR=14; PT=36; PB=42
  ww <- W-PL-PR; h <- H-PT-PB
  r  <- if (arms == 2) k / (1 + k) else 1
  f  <- 1 - dr / 100
  ts_vec <- seq(0, L, length.out=161)
  pts <- sapply(ts_vec, function(t) {
    eC <- if (arms == 2) N * (1-r) * f * cumP_fn(t, l0, A) else 0
    N * r * f * cumP_fn(t, l1, A) + eC
  })
  eMax <- pts[length(pts)]
  yMax <- max(1, eMax) * 1.14
  xS   <- function(t) PL + (t / L) * ww
  yS   <- function(e) PT + h - (e / yMax) * h
  pth  <- paste0("M", paste0(sprintf("%.1f", xS(ts_vec)), ",", sprintf("%.1f", yS(pts)), collapse="L"))
  
  yStep <- max(1, ceiling(yMax / 4 / 5) * 5)
  yTk   <- Filter(function(y) y <= yMax * 1.05, (0:4) * yStep)
  xTk   <- round(seq(0, 5) / 5 * L)
  
  s <- sprintf('<svg viewBox="0 0 %d %d" style="width:100%%;height:auto">', W, H)
  for (y in yTk)
    s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>',
                           PL, PL+ww, yS(y), yS(y), BDR))
  if (A > 0 && A < L) {
    aw <- xS(A) - xS(0)
    s <- paste0(s, sprintf('<rect x="%.1f" y="%d" width="%.1f" height="%d" fill="%s10"/>',
                           xS(0), PT, aw, h, PK))
    s <- paste0(s, sprintf('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="1" stroke-dasharray="4,3"/>',
                           xS(A), xS(A), PT, PT+h, PK))
    s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="8" fill="%s">A=%g %s</text>',
                           xS(A), PT+11, PK, A, unit_short))
  }
  if (eMax > 0 && eMax / yMax <= 1.02) {
    s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="%s" stroke-width="1" stroke-dasharray="4,3"/>',
                           PL, PL+ww, yS(eMax), yS(eMax), GN))
    s <- paste0(s, sprintf('<text x="%.1f" y="%.1f" text-anchor="end" font-size="8" fill="%s">D~%g</text>',
                           PL+ww-4, yS(eMax)-4, GN, round(eMax)))
  }
  s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="%s" stroke-width="1.5"/>',
                         PL, PL+ww, PT+h, PT+h, GY))
  s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="%s" stroke-width="1.5"/>',
                         PL, PL, PT, PT+h, GY))
  s <- paste0(s, sprintf('<path d="%s" stroke="%s" stroke-width="2" fill="none"/>', pth, GY))
  for (x in xTk)
    s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="9" fill="%s">%g</text>',
                           xS(x), PT+h+14, MU, x))
  for (y in yTk)
    s <- paste0(s, sprintf('<text x="%d" y="%.1f" text-anchor="end" font-size="9" fill="%s">%g</text>',
                           PL-5, yS(y)+4, MU, y))
  s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="9" fill="%s">Calendar time (%s), L=%g</text>',
                         PL+ww/2, H-3, MU, unit_long, L))
  s <- paste0(s, sprintf('<text x="12" y="%.1f" text-anchor="middle" font-size="9" fill="%s" transform="rotate(-90,12,%.1f)">Cum. events</text>',
                         PT+h/2, MU, PT+h/2))
  paste0(s, "</svg>")
}

make_heat_svg <- function(l0, A, L, dr, al, sd, cHR, cN, arms, k) {
  HRv <- c(0.40,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95)
  Nv  <- if (arms == 1) c(10,20,30,50,75,100,150,200,300,400,500) else
    c(20,40,60,80,100,150,200,250,300,400,500)
  W=550; H=300; PL=50; PR=80; PT=34; PB=46
  ww <- W-PL-PR; hh <- H-PT-PB
  cW <- ww / length(Nv); cH <- hh / length(HRv)
  xS <- function(ci) PL + ci * cW
  yS <- function(ri) PT + ri * cH
  pc <- function(p) {
    p <- max(0, min(1, p))
    if (p < 0.5) {
      t_ <- p / 0.5
      sprintf("rgb(%d,%d,%d)",
              round(225-30*t_), round(75+145*t_), round(75+145*t_))
    } else {
      t_ <- (p - 0.5) / 0.5
      sprintf("rgb(%d,%d,%d)",
              round(195-170*t_), round(220-60*t_), round(220+35*t_))
    }
  }
  cNi <- which.min(abs(Nv - cN)) - 1  # 0-based
  cHi <- which.min(abs(HRv - cHR)) - 1
  
  s <- sprintf('<svg viewBox="0 0 %d %d" style="width:100%%;height:auto">', W, H)
  
  for (ri in seq_along(HRv) - 1) {
    hr <- HRv[ri+1]
    for (ci in seq_along(Nv) - 1) {
      n   <- Nv[ci+1]
      l1  <- l0 * hr
      r   <- if (arms == 2) k / (1 + k) else NULL
      pe  <- if (arms == 1)
        pEv_fn(l1, A, L) * (1 - dr/100)
      else
        (r * pEv_fn(l1, A, L) + (1-r) * pEv_fn(l0, A, L)) * (1 - dr/100)
      p   <- pwrN_fn(n, hr, al, pe, sd, r)
      hi  <- (ri == cHi && ci == cNi)
      s <- paste0(s, sprintf(
        '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" stroke="none"/>',
        xS(ci), yS(ri), cW, cH, pc(p)))
      
      if (hi) {
        s <- paste0(s, sprintf(
          '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="4" ry="4" fill="none" stroke="%s" stroke-width="3"/>',
          xS(ci)+1.5, yS(ri)+1.5, cW-3, cH-3, BG))
        s <- paste0(s, sprintf(
          '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="4" ry="4" fill="none" stroke="%s" stroke-width="1.8"/>',
          xS(ci)+1.5, yS(ri)+1.5, cW-3, cH-3, "#ffffff"))
      }
      tx <- xS(ci) + cW/2; ty <- yS(ri) + cH/2 + 4
      s <- paste0(s, sprintf(
        '<text x="%.1f" y="%.1f" text-anchor="middle" font-size="%d" fill="%s" font-weight="%d">%d%%</text>',
        tx, ty, if(cW > 40) 9 else 7,
        if(p > 0.55) "#0f172a" else "#1e3a5f",
        if(hi) 700 else 400, round(p*100)))
    }
  }
  s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="%s" stroke-width="1.5"/>',
                         PL, PL+ww, PT+hh, PT+hh, GY))
  s <- paste0(s, sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="%s" stroke-width="1.5"/>',
                         PL, PL, PT, PT+hh, GY))
  for (ci in seq_along(Nv) - 1)
    s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="8" fill="%s">%d</text>',
                           xS(ci)+cW/2, PT+hh+14, MU, Nv[ci+1]))
  for (ri in seq_along(HRv) - 1)
    s <- paste0(s, sprintf('<text x="%d" y="%.1f" text-anchor="end" font-size="8" fill="%s">%.2f</text>',
                           PL-4, yS(ri)+cH/2+4, MU, HRv[ri+1]))
  s <- paste0(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="9" fill="%s">Sample Size N</text>',
                         PL+ww/2, H-4, MU))
  s <- paste0(s, sprintf('<text x="12" y="%.1f" text-anchor="middle" font-size="9" fill="%s" transform="rotate(-90,12,%.1f)">HR</text>',
                         PT+hh/2, MU, PT+hh/2))
  for (i in 0:39) {
    s <- paste0(s, sprintf('<rect x="%d" y="%.1f" width="14" height="%.1f" fill="%s"/>',
                           W-PR+10, PT+i*(hh/40), hh/40+1, pc(1-i/39)))
  }
  s <- paste0(s, sprintf('<text x="%d" y="%d" text-anchor="middle" font-size="7" fill="%s">100%%</text>',
                         W-PR+17, PT-3, MU))
  s <- paste0(s, sprintf('<text x="%d" y="%d" text-anchor="middle" font-size="7" fill="%s">0%%</text>',
                         W-PR+17, PT+hh+12, MU))
  paste0(s, "</svg>")
}

# ─── CSS helpers ──────────────────────────────────────────────────────────────
css <- function() {
  tags$style(HTML(paste0("
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
    html, body { height: 100%; background: ", BG, "; color: ", TX, ";
      font-family: 'Segoe UI', system-ui, sans-serif; }
    ::-webkit-scrollbar { width: 4px; }
    ::-webkit-scrollbar-track { background: ", PAN, "; }
    ::-webkit-scrollbar-thumb { background: ", B2, "; border-radius: 2px; }
    input[type=range] { height: 3px; accent-color: ", BL, "; width: 100%; }
    input[type=number], input[type=text] {
      background: ", BG, "; border: 1px solid ", DM, ";
      border-radius: 6px; color: ", TX, "; padding: 5px 9px;
      font-size: 13px; font-family: monospace; outline: none; width: 100%; }
    input[type=number]:focus, input[type=text]:focus { border-color: ", BL, "; }
    .pill-btn {
      padding: 4px 11px; border-radius: 20px; cursor: pointer;
      font-family: inherit; font-size: 12px; border: 1px solid ", DM, ";
      background: ", PAN, "; color: ", MU, "; margin-right: 4px; margin-bottom: 4px;
      display: inline-block; }
    .pill-btn.active-BL { border-color: ", BL, "; background: ", BL, "22; color: ", BL, "; }
    .pill-btn.active-TL { border-color: ", TL, "; background: ", TL, "22; color: ", TL, "; }
    .pill-btn.active-GN { border-color: ", GN, "; background: ", GN, "22; color: ", GN, "; }
    .pill-btn.active-OR { border-color: ", OR, "; background: ", OR, "22; color: ", OR, "; }
    .pill-btn.active-GY { border-color: ", GY, "; background: ", GY, "22; color: ", GY, "; }
    .sec-hdr {
      font-size: 10px; letter-spacing: 2px; color: ", DM, ";
      text-transform: uppercase; border-bottom: 1px solid ", BDR, ";
      padding-bottom: 5px; margin-bottom: 10px; margin-top: 12px; }
    .kpi-card {
      background: ", BG, "; border-radius: 8px; padding: 10px 13px;
      flex: 1 1 90px; min-width: 80px; }
    .kpi-label { color: ", MU, "; font-size: 10px; letter-spacing: 1px;
      text-transform: uppercase; margin-bottom: 4px; }
    .kpi-val { font-size: 20px; font-weight: 700; font-family: monospace; line-height: 1; }
    .kpi-sub { color: ", MU, "; font-size: 11px; margin-top: 3px; }
    .viz-card {
      background: ", PAN, "; border: 1px solid ", BDR, ";
      border-radius: 10px; padding: 13px 16px; margin-bottom: 14px; }
    .viz-card-hdr {
      display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px; }
    .toggle-btn {
      background: none; border: 1px solid ", DM, "; border-radius: 5px;
      color: ", MU, "; cursor: pointer; font-size: 11px; padding: 2px 8px; }
    .panel-box {
      background: ", PAN, "; border: 1px solid ", BDR, ";
      border-radius: 10px; padding: 13px 16px; margin-bottom: 14px; }
    table.sens-table { width: 100%; border-collapse: collapse; font-size: 12px; }
    table.sens-table th {
      padding: 5px 9px; color: ", MU, "; text-align: left;
      border-bottom: 1px solid ", BDR, "; font-weight: 500; font-size: 11px;
      white-space: nowrap; }
    table.sens-table td { padding: 5px 9px; font-family: monospace; }
    .nav-tabs > li > a { background: ", PAN, " !important; color: ", MU, " !important;
      border-color: ", BDR, " !important; }
    .nav-tabs > li.active > a, .nav-tabs > li.active > a:focus,
    .nav-tabs > li.active > a:hover {
      background: ", B2, " !important; color: ", TX, " !important;
      border-color: ", B2, " !important; }
    .tab-content { background: ", BG, "; }
    .sidebar-panel { width: 298px !important; background: #05101e !important;
      border-right: 1px solid ", BDR, "; padding: 0 0 32px; min-height: 100vh; }
    .main-panel { background: ", BG, "; }
    ol li { margin-bottom: 3px; }
    a { color: ", BL, "; }
  ")))
}

# ─── UI ───────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  useShinyjs(),
  css(),
  title = "TTE Sample Size Calculator",
  div(style=paste0("display:flex; min-height:100vh; background:", BG, ";"),
      
      # ── Sidebar ──────────────────────────────────────────────────────────────
      div(class="sidebar-panel", style="overflow-y:auto; flex-shrink:0;",
          div(style=paste0("padding:13px 15px; border-bottom:1px solid ", BDR, "; margin-bottom:12px;"),
              div(style="display:flex; align-items:center; gap:8px; margin-bottom:2px;",
                  div(style=paste0("width:7px; height:7px; border-radius:50%; background:", BL,
                                   "; box-shadow:0 0 8px ", BL, "88;")),
                  span(style="font-weight:700; font-size:13px;", "Time-to-Event Sample Size")
              ),
              div(style=paste0("color:", DM, "; font-size:11px;"),
                  "Exponential \u00b7 Log-rank \u00b7 Schoenfeld / Lachin-Foulkes")
          ),
          div(style="padding: 0 14px;",
              div(class="sec-hdr", "Number of Arms"),
              div(
                actionButton("btn_arms1", "Single Arm", class="pill-btn"),
                actionButton("btn_arms2", "Two Arms",   class="pill-btn")
              ),
              
              # Two-arm extras (conditionally shown via shinyjs)
              hidden(div(id="two_arm_extras",
                         div(class="sec-hdr", "Test Objective"),
                         div(
                           actionButton("btn_obj_eq",  "Equality",        class="pill-btn"),
                           actionButton("btn_obj_ni",  "Non-inferiority", class="pill-btn"),
                           actionButton("btn_obj_sup", "Superiority",     class="pill-btn"),
                           actionButton("btn_obj_eqv", "Equivalence",     class="pill-btn")
                         ),
                         hidden(div(id="dlt_box",
                                    div(style=paste0("color:", GY, "; font-size:12px; margin-bottom:3px;"),
                                        textOutput("dlt_lbl", inline=TRUE)),
                                    numericInput("dlt", NULL, value=0.2, min=0.001, max=2, step=0.01)
                         )),
                         div(
                           style = paste0("color:", GY, "; font-size:12px; margin-bottom:3px;"),
                           HTML("Allocation ratio k = n<sub>t</sub> / n<sub>c</sub>")
                         ),
                         numericInput("k", NULL, value = 1, min = 0.1, max = 5, step = 0.1),
                         div(
                           style = paste0("color:", MU, "; font-size:10px; font-family:monospace; margin-top:-6px; margin-bottom:8px;"),
                           textOutput("k_note", inline = TRUE)
                         )
              )),
              
              div(class="sec-hdr", "Test Direction"),
              div(
                actionButton("btn_sd1", "1-Sided", class="pill-btn"),
                actionButton("btn_sd2", "2-Sided", class="pill-btn")
              ),
              
              div(class="sec-hdr", "Primary Endpoint"),
              div(
                actionButton("btn_ep_PFS",   "PFS",   class="pill-btn"),
                actionButton("btn_ep_OS",    "OS",    class="pill-btn"),
                actionButton("btn_ep_DFS",   "DFS",   class="pill-btn"),
                actionButton("btn_ep_EFS",   "EFS",   class="pill-btn"),
                actionButton("btn_ep_Other", "Other", class="pill-btn")
              ),
              hidden(div(id="ep_other_box",
                         textInput("ep_other", NULL, placeholder="Endpoint name\u2026"),
                         style="margin-bottom:8px;"
              )),
              
              div(class="sec-hdr", "Time Unit"),
              div(
                actionButton("btn_tu_mo", "Months", class="pill-btn"),
                actionButton("btn_tu_yr", "Years",  class="pill-btn")
              ),
              
              div(class="sec-hdr", "Input Format"),
              div(
                actionButton("btn_im_med",  "Median",         class="pill-btn"),
                actionButton("btn_im_St",   "Landmark S(t*)", class="pill-btn"),
                actionButton("btn_im_haz",  "Hazard rate",    class="pill-btn")
              ),
              hidden(div(id="ts_box",
                         div(style=paste0("color:", GY, "; font-size:12px; margin-bottom:3px;"),
                             textOutput("ts_lbl", inline=TRUE)),
                         numericInput("ts", NULL, value=18, min=1, max=120, step=0.5)
              )),
              
              div(class="sec-hdr", textOutput("hyp_lbl", inline=TRUE)),
              
              # Median inputs
              div(id="im_median_box",
                  div(style=paste0("color:", GY, "; font-size:12px; margin-bottom:3px;"),
                      textOutput("mc_lbl", inline=TRUE)),
                  numericInput("mc", NULL, value=18, min=0.01, max=240, step=0.5),
                  div(style=paste0("color:", MU, "; font-size:10px; font-family:monospace; margin-top:-6px; margin-bottom:6px;"),
                      textOutput("note0", inline=TRUE)),
                  div(style=paste0("color:", GY, "; font-size:12px; margin-bottom:3px;"),
                      textOutput("mt_lbl", inline=TRUE)),
                  numericInput("mt", NULL, value=27, min=0.01, max=240, step=0.5),
                  div(style=paste0("color:", MU, "; font-size:10px; font-family:monospace; margin-top:-6px; margin-bottom:6px;"),
                      textOutput("note1", inline=TRUE))
              ),
              # Landmark S(t*) inputs
              hidden(div(id="im_St_box",
                         div(style=paste0("color:", GY, "; font-size:12px; margin-bottom:3px;"),
                             textOutput("s0_lbl", inline=TRUE)),
                         numericInput("s0", NULL, value=0.46, min=0.01, max=0.99, step=0.01),
                         div(style=paste0("color:", MU, "; font-size:10px; font-family:monospace; margin-top:-6px; margin-bottom:6px;"),
                             textOutput("note0_St", inline=TRUE)),
                         div(style=paste0("color:", GY, "; font-size:12px; margin-bottom:3px;"),
                             textOutput("s1_lbl", inline=TRUE)),
                         numericInput("s1", NULL, value=0.63, min=0.01, max=0.99, step=0.01),
                         div(style=paste0("color:", MU, "; font-size:10px; font-family:monospace; margin-top:-6px; margin-bottom:6px;"),
                             textOutput("note1_St", inline=TRUE))
              )),
              # Hazard rate inputs
              hidden(div(id="im_hazard_box",
                         div(style=paste0("color:", GY, "; font-size:12px; margin-bottom:3px;"),
                             textOutput("lc_lbl", inline=TRUE)),
                         numericInput("lc", NULL, value=0.0385, min=0.0001, max=5, step=0.0001),
                         div(style=paste0("color:", MU, "; font-size:10px; font-family:monospace; margin-top:-6px; margin-bottom:6px;"),
                             textOutput("note0_haz", inline=TRUE)),
                         div(style=paste0("color:", GY, "; font-size:12px; margin-bottom:3px;"),
                             textOutput("lt_lbl", inline=TRUE)),
                         numericInput("lt", NULL, value=0.0257, min=0.0001, max=5, step=0.0001),
                         div(style=paste0("color:", MU, "; font-size:10px; font-family:monospace; margin-top:-6px; margin-bottom:6px;"),
                             textOutput("note1_haz", inline=TRUE))
              )),
              
              # HR display box
              div(style=paste0("background:", BG, "; border-radius:7px; padding:8px 10px;",
                               "margin-bottom:12px; font-size:11px; font-family:monospace;",
                               "line-height:2; border: 1px solid ", DM, ";"),
                  htmlOutput("hr_box")),
              
              div(class="sec-hdr", "Type I Error & Power"),
              div(style=paste0("display:flex; justify-content:space-between; margin-bottom:3px;"),
                  span(style=paste0("color:", GY, "; font-size:12px;"),
                       textOutput("al_lbl", inline=TRUE)),
                  span(style=paste0("color:", TX, "; font-size:12px; font-family:monospace;"),
                       textOutput("al_val", inline=TRUE))
              ),
              sliderInput("al", NULL, min=0.01, max=0.20, step=0.01, value=0.05, ticks=FALSE),
              div(style=paste0("display:flex; justify-content:space-between; margin-bottom:3px;"),
                  span(style=paste0("color:", GY, "; font-size:12px;"), "Power (1 \u2212 \u03b2)"),
                  span(style=paste0("color:", TX, "; font-size:12px; font-family:monospace;"),
                       textOutput("pw_val", inline=TRUE))
              ),
              sliderInput("pw", NULL, min=0.70, max=0.99, step=0.01, value=0.80, ticks=FALSE),
              
              div(class="sec-hdr", "Timing"),
              div(style=paste0("display:flex; justify-content:space-between; margin-bottom:3px;"),
                  span(style=paste0("color:", GY, "; font-size:12px;"),
                       textOutput("A_lbl", inline=TRUE)),
                  span(style=paste0("color:", TX, "; font-size:12px; font-family:monospace;"),
                       textOutput("A_val", inline=TRUE))
              ),
              sliderInput("A", NULL, min=0.5, max=60, step=0.5, value=24, ticks=FALSE),
              div(style=paste0("display:flex; justify-content:space-between; margin-bottom:3px;"),
                  span(style=paste0("color:", GY, "; font-size:12px;"),
                       textOutput("L_lbl", inline=TRUE)),
                  span(style=paste0("color:", TX, "; font-size:12px; font-family:monospace;"),
                       textOutput("L_val", inline=TRUE))
              ),
              sliderInput("L", NULL, min=1.0, max=120, step=0.5, value=36, ticks=FALSE),
              div(style=paste0("font-size:11px; color:", DM, "; margin-top:-5px; margin-bottom:10px;"),
                  textOutput("F_note", inline=TRUE)),
              
              div(class="sec-hdr", "Dropout / Loss to Follow-up"),
              div(style=paste0("display:flex; justify-content:space-between; margin-bottom:3px;"),
                  span(style=paste0("color:", GY, "; font-size:12px;"), "Anticipated dropout (%)"),
                  span(style=paste0("color:", TX, "; font-size:12px; font-family:monospace;"),
                       textOutput("dr_val", inline=TRUE))
              ),
              sliderInput("dr", NULL, min=0, max=40, step=1, value=10, ticks=FALSE),
              div(style=paste0("background:", BG, "; border:1px solid ", BDR, "; border-radius:6px;",
                               "padding:7px 10px; font-size:11px; color:", MU, "; line-height:1.7;"),
                  HTML(paste0("The <strong>total proportion</strong> of enrolled patients not contributing ",
                              "an observed event due to withdrawal or loss to follow-up, ",
                              "<strong>not an annual rate</strong>.<br>",
                              "<span style='font-family:monospace; font-size:10px; color:", GY, ";'>",
                              "p\u2091 = p(event) \u00d7 (1 \u2212 dropout)</span><br>",
                              "<span style='font-size:10px; color:", MU, ";'>",
                              "p\u2091 = dropout-adjusted event probability used to derive N from D.</span>"))
              )
          )
      ),
      
      # ── Main panel ────────────────────────────────────────────────────────────
      div(class="main-panel", style="flex:1; overflow-y:auto; padding:18px 22px;",
          div(style=paste0("color:", DM, "; font-size:11px; letter-spacing:2px;",
                           "text-transform:uppercase; margin-bottom:12px;"),
              textOutput("subtitle")),
          
          # L <= A warning
          div(id="la_warn",
              div(style=paste0("background:#2d0a0a; border:1px solid ", RD, ";",
                               "border-radius:8px; padding:10px 14px; margin-bottom:14px;",
                               "color:", RD, "; font-size:13px;"),
                  textOutput("la_warn_txt"))
          ),
          
          # KPI row
          div(style="display:flex; gap:8px; flex-wrap:wrap; margin-bottom:16px;",
              htmlOutput("kpi_row")),
          
          # Interpretation box
          htmlOutput("interp_box"),
          
          # Step-by-step calculation
          htmlOutput("calc_steps"),
          
          # ── Visualizations in tabs ────────────────────────────────────────────
          div(class="panel-box",
              tabsetPanel(id="viz_tabs",
                          tabPanel("1 \u00b7 Survival Curves",
                                   div(style="margin-top:12px;",
                                       htmlOutput("svg_surv"),
                                       div(style=paste0("font-size:11px; color:", MU, "; margin-top:7px; line-height:1.6;"),
                                           textOutput("cap_surv"))
                                   )
                          ),
                          tabPanel("2 \u00b7 Test-Statistic Distribution",
                                   div(style="margin-top:12px;",
                                       htmlOutput("svg_z"),
                                       div(style=paste0("font-size:11px; color:", MU, "; margin-top:7px; line-height:1.6;"),
                                           textOutput("cap_z"))
                                   )
                          ),
                          tabPanel("3 \u00b7 Cumulative Events",
                                   div(style="margin-top:12px;",
                                       htmlOutput("svg_evt"),
                                       div(style=paste0("font-size:11px; color:", MU, "; margin-top:7px; line-height:1.6;"),
                                           textOutput("cap_evt"))
                                   )
                          ),
                          tabPanel("4 \u00b7 Power Heatmap",
                                   div(style="margin-top:12px;",
                                       htmlOutput("svg_heat"),
                                       div(style=paste0("font-size:11px; color:", MU, "; margin-top:7px; line-height:1.6;"),
                                           textOutput("cap_heat"))
                                   )
                          )
              )
          ),
          
          # Sensitivity table
          div(class="panel-box",
              div(class="sec-hdr", "Sensitivity: Events & Sample Size Around Your Design"),
              div(style=paste0("font-size:11px; color:", MU, "; margin-bottom:10px; line-height:1.8;"),
                  htmlOutput("sens_note")),
              div(style="overflow-x:auto;",
                  uiOutput("sens_table"))
          ),
          
          # Budget back-calculation
          div(class="panel-box",
              div(class="sec-hdr", "Budget Constraint: Back-Calculate from N"),
              div(style=paste0("font-size:12px; color:", MU, "; margin-bottom:12px; line-height:1.7;"),
                  "Enter your maximum feasible sample size and we will calculate what ",
                  "power and \u03b1 are achievable given your H\u2080 / H\u2081 and study design."),
              div(style="display:flex; gap:12px; align-items:center; margin-bottom:14px;",
                  div(style="flex:0 0 200px;",
                      div(style=paste0("color:", GY, "; font-size:12px; margin-bottom:4px;"),
                          textOutput("budget_lbl", inline=TRUE)),
                      numericInput("budget_n", NULL, value=NA, min=1, max=5000, step=1)
                  )
              ),
              htmlOutput("budget_result")
          ),
          
          # References
          div(class="panel-box",
              div(class="sec-hdr", "References"),
              htmlOutput("references")
          )
      )
  )
)

# ─── Server ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── Reactive state (mirrors JSX useState hooks) ───────────────────────────
  rv <- reactiveValues(
    arms = 1,
    obj  = "equality",
    sd   = 1,
    ep   = "PFS",
    tu   = "months",
    im   = "median"
  )
  
  # ── Pill button handlers ──────────────────────────────────────────────────
  observeEvent(input$btn_arms1, { rv$arms <- 1; toggle_arm_ui(1) })
  observeEvent(input$btn_arms2, { rv$arms <- 2; toggle_arm_ui(2) })
  
  toggle_arm_ui <- function(a) {
    if (a == 2) shinyjs::show("two_arm_extras") else shinyjs::hide("two_arm_extras")
    update_pill_classes()
  }
  
  update_pill_classes <- function() {
    #Remove all active classes from all pill buttons
    ids <- c(
      "btn_arms1","btn_arms2",
      "btn_obj_eq","btn_obj_ni","btn_obj_sup","btn_obj_eqv",
      "btn_sd1","btn_sd2",
      "btn_ep_PFS","btn_ep_OS","btn_ep_DFS","btn_ep_EFS","btn_ep_Other",
      "btn_tu_mo","btn_tu_yr",
      "btn_im_med","btn_im_St","btn_im_haz"
    )
    
    active_classes <- c("active-BL","active-TL","active-GN","active-OR","active-GY")
    
    for (id in ids) {
      for (cls in active_classes) {
        shinyjs::removeClass(selector = paste0("#", id), class = cls)
      }
    }
    
    #Apply active class to current selections
    shinyjs::addClass(
      id = if (rv$arms == 1) "btn_arms1" else "btn_arms2",
      class = "active-BL"
    )
    
    shinyjs::addClass(
      id = switch(rv$obj,
                  "equality"       = "btn_obj_eq",
                  "noninferiority" = "btn_obj_ni",
                  "superiority"    = "btn_obj_sup",
                  "equivalence"    = "btn_obj_eqv"
      ),
      class = "active-TL"
    )
    
    shinyjs::addClass(
      id = if (rv$sd == 1) "btn_sd1" else "btn_sd2",
      class = "active-GN"
    )
    
    shinyjs::addClass(
      id = switch(rv$ep,
                  "PFS"   = "btn_ep_PFS",
                  "OS"    = "btn_ep_OS",
                  "DFS"   = "btn_ep_DFS",
                  "EFS"   = "btn_ep_EFS",
                  "Other" = "btn_ep_Other"
      ),
      class = "active-OR"
    )
    
    shinyjs::addClass(
      id = if (rv$tu == "months") "btn_tu_mo" else "btn_tu_yr",
      class = "active-GY"
    )
    
    shinyjs::addClass(
      id = switch(rv$im,
                  "median" = "btn_im_med",
                  "St"     = "btn_im_St",
                  "hazard" = "btn_im_haz"
      ),
      class = "active-BL"
    )
  }
  
  observeEvent(input$btn_obj_eq,  { rv$obj <- "equality";       toggle_obj_ui(); update_pill_classes() })
  observeEvent(input$btn_obj_ni,  { rv$obj <- "noninferiority"; toggle_obj_ui(); update_pill_classes() })
  observeEvent(input$btn_obj_sup, { rv$obj <- "superiority";    toggle_obj_ui(); update_pill_classes() })
  observeEvent(input$btn_obj_eqv, { rv$obj <- "equivalence";    toggle_obj_ui(); update_pill_classes() })
  toggle_obj_ui <- function() {
    if (rv$obj %in% c("noninferiority","superiority","equivalence"))
      shinyjs::show("dlt_box") else shinyjs::hide("dlt_box")
  }
  
  observeEvent(input$btn_sd1, { rv$sd <- 1; update_pill_classes() })
  observeEvent(input$btn_sd2, { rv$sd <- 2; update_pill_classes() })
  
  observeEvent(input$btn_ep_PFS,   { rv$ep <- "PFS"; update_pill_classes() })
  observeEvent(input$btn_ep_OS,    { rv$ep <- "OS";  update_pill_classes() })
  observeEvent(input$btn_ep_DFS,   { rv$ep <- "DFS"; update_pill_classes() })
  observeEvent(input$btn_ep_EFS,   { rv$ep <- "EFS"; update_pill_classes() })
  observeEvent(input$btn_ep_Other, { rv$ep <- "Other"; shinyjs::show("ep_other_box"); update_pill_classes() })
  observe({
    if (rv$ep != "Other") shinyjs::hide("ep_other_box")
  })
  
  observeEvent(input$btn_tu_mo, { rv$tu <- "months"; update_timing_sliders(); update_pill_classes() })
  observeEvent(input$btn_tu_yr, { rv$tu <- "years";  update_timing_sliders(); update_pill_classes() })
  update_timing_sliders <- function() {
    if (rv$tu == "years") {
      updateSliderInput(session,"A", max=10,  step=0.25)
      updateSliderInput(session,"L", max=20,  step=0.25)
    } else {
      updateSliderInput(session,"A", max=60,  step=0.5)
      updateSliderInput(session,"L", max=120, step=0.5)
    }
  }
  
  observeEvent(input$btn_im_med, { rv$im <- "median"; toggle_im_ui(); update_pill_classes() })
  observeEvent(input$btn_im_St,  { rv$im <- "St";     toggle_im_ui(); update_pill_classes() })
  observeEvent(input$btn_im_haz, { rv$im <- "hazard"; toggle_im_ui(); update_pill_classes() })
  toggle_im_ui <- function() {
    shinyjs::toggle("im_median_box", condition=(rv$im == "median"))
    shinyjs::toggle("im_St_box",     condition=(rv$im == "St"))
    shinyjs::toggle("im_hazard_box", condition=(rv$im == "hazard"))
    shinyjs::toggle("ts_box",        condition=(rv$im == "St"))
  }
  
  observe({
    toggle_arm_ui(rv$arms)
    toggle_obj_ui()
    toggle_im_ui()
    update_pill_classes()
  })
  
  # L must be > A — keep slider min in sync
  observeEvent(input$A, {
    minL <- input$A + if (rv$tu == "years") 0.25 else 0.5
    if (!is.null(input$L) && input$L <= input$A)
      updateSliderInput(session, "L", value=minL)
    updateSliderInput(session, "L", min=minL)
  })
  
  # ── Derived unit labels ──────────────────────────────────────────────────
  U   <- reactive(if (rv$tu == "years") "yr" else "mo")
  Ul  <- reactive(if (rv$tu == "years") "years" else "months")
  lU  <- reactive(if (rv$tu == "years") "/yr" else "/mo")
  epN <- reactive(if (rv$ep == "Other") (if(nchar(input$ep_other)>0) input$ep_other else "Endpoint") else rv$ep)
  
  # ── Derive hazard rates from input mode ──────────────────────────────────
  l0 <- reactive({
    switch(rv$im,
           "median" = lm_fn(input$mc),
           "hazard" = input$lc,
           "St"     = ls_fn(input$s0, input$ts)
    )
  })
  l1 <- reactive({
    switch(rv$im,
           "median" = lm_fn(input$mt),
           "hazard" = input$lt,
           "St"     = ls_fn(input$s1, input$ts)
    )
  })
  m0 <- reactive(if (is.finite(l0())) ml_fn(l0()) else NaN)
  m1 <- reactive(if (is.finite(l1())) ml_fn(l1()) else NaN)
  HR <- reactive({
    if (is.finite(l0()) && is.finite(l1()) && l0() > 0) l1()/l0() else NaN
  })
  F_ <- reactive(max(0, input$L - input$A))
  
  # ── Core calculation ──────────────────────────────────────────────────────
  res <- reactive({
    req(is.finite(l0()), is.finite(l1()), l0()>0, l1()>0, input$L > input$A)
    if (rv$arms == 1)
      calcOne(l0(), l1(), input$al, input$pw, rv$sd, input$A, input$L, input$dr)
    else
      calcTwo(l0(), l1(), input$al, input$pw, rv$sd, rv$obj,
              input$dlt, input$k, input$A, input$L, input$dr)
  })
  
  v   <- reactive(!is.null(res()) && is.finite(res()$D) && res()$D > 0)
  cN  <- reactive(if (v()) (if(rv$arms==1) res()$N else res()$Nt) else 100)
  
  # ── Sidebar labels ────────────────────────────────────────────────────────
  output$dlt_lbl  <- renderText(paste0("Margin \u03b4 (log scale): exp(\u03b4) = ",
                                       sprintf("%.3f", exp(input$dlt))))
  output$k_note   <- renderText(paste0("r=", sprintf("%.0f%%", input$k/(1+input$k)*100), " treatment"))
  output$ts_lbl   <- renderText(paste0("Landmark time t* (", Ul(), ")"))
  output$hyp_lbl  <- renderText(if(rv$arms==1) "Null & Alternative Hypotheses" else "Control & Treatment")
  output$mc_lbl   <- renderText(if(rv$arms==1) paste0("H\u2080 benchmark median (", Ul(), ")")
                                else paste0("Control median m\u1D9C (", Ul(), ")"))
  output$mt_lbl   <- renderText(if(rv$arms==1) paste0("H\u2081 target median (", Ul(), ")")
                                else paste0("Treatment median m\u209c (", Ul(), ")"))
  output$s0_lbl   <- renderText(if(rv$arms==1) paste0("H\u2080 survival at t* = ", input$ts, " ", U())
                                else paste0("Control S(t* = ", input$ts, " ", U(), ")"))
  output$s1_lbl   <- renderText(if(rv$arms==1) paste0("H\u2081 survival at t* = ", input$ts, " ", U())
                                else paste0("Treatment S(t* = ", input$ts, " ", U(), ")"))
  output$lc_lbl   <- renderText(if(rv$arms==1) paste0("H\u2080 hazard rate \u03bb\u2080 (", lU(), ")")
                                else paste0("Control hazard rate \u03bb\u1D9C (", lU(), ")"))
  output$lt_lbl   <- renderText(if(rv$arms==1) paste0("H\u2081 hazard rate \u03bb\u2081 (", lU(), ")")
                                else paste0("Treatment hazard rate \u03bb\u209c (", lU(), ")"))
  output$note0    <- renderText({
    if (rv$im == "median") paste0("\u03bb\u2080 = ln(2)/", input$mc, " = ",
                                  if(is.finite(l0())) sprintf("%.4f", l0()) else "\u2014", " ", lU())
    else ""
  })
  output$note1    <- renderText({
    if (rv$im == "median") paste0("\u03bb\u2081 = ln(2)/", input$mt, " = ",
                                  if(is.finite(l1())) sprintf("%.4f", l1()) else "\u2014", " ", lU())
    else ""
  })
  output$note0_St <- renderText({
    if (rv$im == "St") paste0("\u03bb\u2080 = \u2212ln(", input$s0, ")/", input$ts, " = ",
                              if(is.finite(l0())) sprintf("%.4f", l0()) else "\u2014", " ", lU())
    else ""
  })
  output$note1_St <- renderText({
    if (rv$im == "St") paste0("\u03bb\u2081 = \u2212ln(", input$s1, ")/", input$ts, " = ",
                              if(is.finite(l1())) sprintf("%.4f", l1()) else "\u2014", " ", lU())
    else ""
  })
  output$note0_haz <- renderText({
    if (rv$im == "hazard") paste0("Implied m\u2080 = ",
                                  if(is.finite(m0())) sprintf("%.2f", m0()) else "\u2014", " ", U())
    else ""
  })
  output$note1_haz <- renderText({
    if (rv$im == "hazard") paste0("Implied m\u2081 = ",
                                  if(is.finite(m1())) sprintf("%.2f", m1()) else "\u2014", " ", U())
    else ""
  })
  
  output$al_lbl <- renderText(paste0("\u03b1: Type I error (",
                                     if(rv$sd==1) "one-sided" else "two-sided", ")"))
  output$al_val <- renderText(sprintf("%.2f", input$al))
  output$pw_val <- renderText(paste0(sprintf("%.0f%%", input$pw*100)))
  output$A_lbl  <- renderText(paste0("Accrual period A (", Ul(), ")"))
  output$A_val  <- renderText(paste0(input$A, " ", U()))
  output$L_lbl  <- renderText(paste0("Maximum follow-up L = A + F (", Ul(), ")"))
  output$L_val  <- renderText(paste0(input$L, " ", U()))
  output$F_note <- renderText(paste0("Additional follow-up after accrual: F = L \u2212 A = ",
                                     sprintf("%.2f", F_()), " ", Ul()))
  output$dr_val <- renderText(paste0(input$dr, "%"))
  output$budget_lbl <- renderText(paste0("Budget N (maximum ",
                                         if(rv$arms==2) "total " else "", "sample size)"))
  
  # ── HR display box ────────────────────────────────────────────────────────
  output$hr_box <- renderUI({
    hr  <- HR(); m0v <- m0(); m1v <- m1()
    hrc <- if (is.finite(hr) && hr < 1) GN else OR
    hrT <- if (is.finite(hr)) (if(hr<1) " : benefit" else " : harm") else ""
    HTML(paste0(
      '<span style="color:', MU, ';">HR = \u03bb\u2081/\u03bb\u2080 = </span>',
      '<span style="color:', hrc, '; font-weight:600;">',
      if(is.finite(hr)) sprintf("%.4f", hr) else "\u2014", '</span>',
      '<span style="color:', MU, ';">', hrT, '</span><br>',
      '<span style="color:', MU, ';">m\u2080 = ',
      if(is.finite(m0v)) sprintf("%.2f", m0v) else "\u2014", ' ', U(),
      ' \u00b7 m\u2081 = ',
      if(is.finite(m1v)) sprintf("%.2f", m1v) else "\u2014", ' ', U(), '</span>'
    ))
  })
  
  # ── Subtitle ──────────────────────────────────────────────────────────────
  output$subtitle <- renderText({
    paste0(epN(), " \u00b7 ", if(rv$arms==1) "Single-arm" else "Two-arm",
           if(rv$arms==2) paste0(" \u00b7 ", rv$obj) else "",
           " \u00b7 ", if(rv$sd==1) "One-sided" else "Two-sided",
           " \u03b1 = ", input$al)
  })
  
  # ── L <= A warning ────────────────────────────────────────────────────────
  observe({
    if (!is.null(input$L) && !is.null(input$A) && input$L <= input$A)
      shinyjs::show("la_warn")
    else
      shinyjs::hide("la_warn")
  })
  output$la_warn_txt <- renderText({
    paste0("L (", input$L, " ", U(), ") must be greater than A (", input$A, " ", U(), ").")
  })
  
  # ── KPI row ───────────────────────────────────────────────────────────────
  kpi_card <- function(lbl, val, sub, col) {
    paste0(
      '<div class="kpi-card" style="border:1px solid ', col, '33;">',
      '<div class="kpi-label">', lbl, '</div>',
      '<div class="kpi-val" style="color:', col, ';">', val, '</div>',
      '<div class="kpi-sub">', sub, '</div></div>'
    )
  }
  output$kpi_row <- renderUI({
    r <- res(); vv <- v()
    hr <- HR(); m0v <- m0(); m1v <- m1()
    hrc <- if (is.finite(hr) && hr < 1) GN else OR
    
    cards <- paste0(
      kpi_card("Required Events (D)",
               if(vv) r$D else "\u2014",
               if(rv$arms==1) paste0("P(event | H\u2081) = ",
                                     if(vv && is.finite(r$pe)) paste0(sprintf("%.1f%%", r$pe*100)) else "\u2014")
               else paste0("P\u0305(event) = ",
                           if(vv && is.finite(r$pAv)) paste0(sprintf("%.1f%%", r$pAv*100)) else "\u2014"),
               BL)
    )
    if (rv$arms == 1) {
      cards <- paste0(cards,
                      kpi_card("Sample Size (N)",
                               if(vv) r$N else "\u2014",
                               paste0("p\u2091 = ", if(vv&&is.finite(r$peff)) sprintf("%.1f%%",r$peff*100) else "\u2014"),
                               GN))
    } else {
      cards <- paste0(cards,
                      kpi_card("Total N",
                               if(vv) r$Nt else "\u2014",
                               paste0("p\u2091 = ", if(vv&&is.finite(r$peff)) sprintf("%.1f%%",r$peff*100) else "\u2014"),
                               GN),
                      kpi_card("n<sub>t</sub>  (treatment)",
                               if(vv) r$nt else "\u2014",
                               paste0("k = ", input$k), TL),
                      kpi_card("n<sub>c</sub> (control)",
                               if(vv) r$nc else "\u2014",
                               paste0("ratio 1:", input$k), OR))
    }
    cards <- paste0(cards,
                    kpi_card("Hazard Ratio (HR)",
                             if(is.finite(hr)) sprintf("%.4f", hr) else "\u2014",
                             paste0("ln(HR) = ", if(is.finite(hr)) sprintf("%.4f", log(hr)) else "\u2014"),
                             hrc),
                    kpi_card(if(rv$arms==1) "H\u2080 Median (m\u2080)" else "Control Median (m\u1D9C)",
                             if(is.finite(m0v)) paste0(sprintf("%.2f", m0v), " ", U()) else "\u2014",
                             paste0("\u03bb\u2080 = ", if(is.finite(l0())) sprintf("%.4f", l0()) else "\u2014", " ", lU()),
                             GY),
                    kpi_card(if(rv$arms==1) "H\u2081 Median (m\u2081)" else "Treatment Median (m\u209c)",
                             if(is.finite(m1v)) paste0(sprintf("%.2f", m1v), " ", U()) else "\u2014",
                             paste0("\u03bb\u2081 = ", if(is.finite(l1())) sprintf("%.4f", l1()) else "\u2014", " ", lU()),
                             GY)
    )
    HTML(paste0('<div style="display:flex; gap:8px; flex-wrap:wrap;">', cards, '</div>'))
  })
  
  # ── Interpretation box ────────────────────────────────────────────────────
  output$interp_box <- renderUI({
    if (!v()) return(NULL)
    r <- res(); hr <- HR()
    hrDir <- if (hr < 1) "reduction" else "increase"
    hrPct <- sprintf("%.0f%%", abs((1-hr)*100))
    effLabel <- switch(rv$im,
                       "median" = paste0("a median ", epN(), " of ", sprintf("%.2f", m1()), " ", Ul(),
                                         " (vs ", sprintf("%.2f", m0()), " ", Ul(), " under H\u2080)"),
                       "hazard" = paste0("a hazard rate of \u03bb\u2081 = ", sprintf("%.4f", l1()), " ", lU(),
                                         " (vs \u03bb\u2080 = ", sprintf("%.4f", l0()), " ", lU(), ")"),
                       "St"     = paste0("a ", input$ts, " ", U(), " ", epN(), " rate of ",
                                         sprintf("%.0f%%", input$s1*100), "% (vs ",
                                         sprintf("%.0f%%", input$s0*100), "% under H\u2080)")
    )
    armDesc <- if (rv$arms==1)
      paste0("This single-arm study compares observed ", epN(), " against a historical benchmark.")
    else
      paste0("This randomised two-arm study compares treatment against control on ", epN(), ".")
    testDesc <- if (rv$arms==1)
      paste0("Using a one-sample log-rank test (",
             if(rv$sd==1) "one" else "two", "-sided \u03b1 = ", input$al, ")")
    else
      paste0("Using a log-rank test (",
             if(rv$sd==1) "one" else "two", "-sided \u03b1 = ", input$al, ", ", rv$obj, ")")
    nDesc <- if (rv$arms==1)
      paste0("A total of ", r$N, " participants are required.")
    else
      paste0("A total of ", r$Nt, " participants are required (",
             r$nt, " treatment, ", r$nc, " control; ratio ", input$k, ":1).")
    
    hrc <- if (hr < 1) GN else OR
    HTML(paste0(
      '<div style="background:#0d1f38; border:1px solid ', B2, ';',
      'border-radius:10px; padding:14px 18px; margin-bottom:16px;',
      'font-size:13px; color:', TX, '; line-height:1.9;">',
      '<span style="color:', GY, '; font-weight:600;">Interpretation: </span>',
      armDesc, ' The design targets ', effLabel, ', corresponding to a ',
      '<span style="color:', hrc, '; font-weight:600;">',
      hrPct, ' ', hrDir, ' in hazard (HR = ', sprintf("%.4f", hr), ')</span>. ',
      testDesc, ', with ', sprintf("%.0f%%", input$pw*100), ' power, ',
      '<span style="color:', BL, '; font-weight:600;">', r$D, ' events</span> ',
      'are required to detect this difference. ',
      'Accounting for uniform accrual over ', input$A, ' ', Ul(),
      if(input$dr > 0) paste0(' and ', input$dr, '% dropout') else '', ', ',
      '<span style="color:', GN, '; font-weight:600;">', nDesc, '</span>',
      '</div>'
    ))
  })
  
  # ── Step-by-step calculation ──────────────────────────────────────────────
  output$calc_steps <- renderUI({
    if (!v()) return(NULL)
    r <- res()
    hr <- HR()
    step1 <- switch(rv$im,
                    "median" = paste0(
                      'Convert medians \u2192 hazard rates: \u03bb = ln(2) / m<br>',
                      '<span style="font-family:monospace; color:', BL, ';">\u03bb\u2080 = ln(2) / ',
                      sprintf("%.2f", m0()), ' = ', sprintf("%.6f", l0()), ' ', lU(), '</span>',
                      ' &nbsp;\u00b7&nbsp; ',
                      '<span style="font-family:monospace; color:', OR, ';">\u03bb\u2081 = ln(2) / ',
                      sprintf("%.2f", m1()), ' = ', sprintf("%.6f", l1()), ' ', lU(), '</span>'),
                    "hazard" = paste0(
                      'Input hazard rates (medians derived: m = ln(2) / \u03bb)<br>',
                      '<span style="font-family:monospace; color:', BL, ';">\u03bb\u2080 = ',
                      sprintf("%.6f", l0()), ' ', lU(), ' \u2192 m\u2080 = ', sprintf("%.2f", m0()), ' ', U(), '</span>',
                      ' &nbsp;\u00b7&nbsp; ',
                      '<span style="font-family:monospace; color:', OR, ';">\u03bb\u2081 = ',
                      sprintf("%.6f", l1()), ' ', lU(), ' \u2192 m\u2081 = ', sprintf("%.2f", m1()), ' ', U(), '</span>'),
                    "St" = paste0(
                      'Convert landmark survivals \u2192 hazard rates: \u03bb = \u2212ln(S) / t*<br>',
                      '<span style="font-family:monospace; color:', BL, ';">\u03bb\u2080 = \u2212ln(',
                      input$s0, ') / ', input$ts, ' = ', sprintf("%.6f", l0()), ' ', lU(), '</span>',
                      ' &nbsp;\u00b7&nbsp; ',
                      '<span style="font-family:monospace; color:', OR, ';">\u03bb\u2081 = \u2212ln(',
                      input$s1, ') / ', input$ts, ' = ', sprintf("%.6f", l1()), ' ', lU(), '</span>')
    )
    step2 <- paste0(
      'Compute HR = \u03bb\u2081 / \u03bb\u2080 = ',
      sprintf("%.4f", l1()), ' / ', sprintf("%.4f", l0()),
      ' = <span style="font-family:monospace;">',
      sprintf("%.4f", hr), '</span><br>',
      '<span style="font-family:monospace; color:', MU, ';">',
      'ln(HR) = ', sprintf("%.4f", log(hr)), '</span>'
    )
    
    effLHR <- if (rv$arms == 1) abs(r$lHR) else abs(if(!is.null(r$eLHR)) r$eLHR else r$lHR)
    r_val  <- if (rv$arms == 2) input$k/(1+input$k) else NULL
    step3 <- if (rv$arms == 1) paste0(
      'Required events: one-sample Wald test on log-hazard scale (',
      if(rv$sd==1) "one" else "two", '-sided \u03b1 = ', input$al, '):<br>',
      '<span style="font-family:monospace; color:', MU, ';">',
      'D = (z\u2081\u208b\u03b1 + z\u2081\u208b\u03b2)\u00b2 / [ln(HR)]\u00b2</span><br>',
      '<span style="font-family:monospace;">D = (',
      sprintf("%.4f", r$z1a), ' + ', sprintf("%.4f", r$z1b), ')\u00b2 / (',
      sprintf("%.4f", r$lHR), ')\u00b2 = <strong style="color:', BL, ';">', r$D, '</strong></span>'
    ) else paste0(
      'Required events: log-rank test (Schoenfeld 1981/1983), r = k/(1+k) = ',
      sprintf("%.4f", r_val), ':<br>',
      '<span style="font-family:monospace; color:', MU, ';">',
      'D = (z\u2081\u208b\u03b1 + z\u2081\u208b\u03b2)\u00b2 / [\u03942 \u00b7 r\u00b7(1\u2212r)]</span><br>',
      '<span style="font-family:monospace;">D = (',
      sprintf("%.4f", r$z1a), ' + ', sprintf("%.4f", r$z1b), ')\u00b2 / (',
      sprintf("%.4f", effLHR), '\u00b2 \u00d7 ',
      sprintf("%.4f", r_val), ' \u00d7 ', sprintf("%.4f", 1-r_val),
      ') = <strong style="color:', BL, ';">', r$D, '</strong></span>'
    )
    
    step4 <- paste0(
      'Per-participant event probability under H\u2081 ',
      '(uniform accrual A = ', input$A, ' ', U(), ', max follow-up L = ', input$L, ' ', U(),
      ', F = ', sprintf("%.2f", F_()), ' ', U(), '):<br>',
      '<span style="font-family:monospace; color:', MU, ';">',
      'p(event) = 1 \u2212 [exp(\u2212\u03bb\u2081\u00b7F) \u2212 exp(\u2212\u03bb\u2081\u00b7L)] / (A\u00b7\u03bb\u2081)</span>'
    )
    if (rv$arms == 1) {
      step4 <- paste0(step4, ' = <span style="font-family:monospace;"><strong>',
                      sprintf("%.4f", r$pe), '</strong></span>')
    } else {
      step4 <- paste0(step4, '<br>',
                      '<span style="font-family:monospace;">p(event | control) = ',
                      sprintf("%.4f", r$pc),
                      ' &nbsp;\u00b7&nbsp; p(event | treatment) = ', sprintf("%.4f", r$pt), '<br>',
                      'Weighted avg: r\u00b7p_t + (1\u2212r)\u00b7p_c = ',
                      sprintf("%.4f", r$pAv), '</span>')
    }
    
    pBase <- if (rv$arms==1) r$pe else r$pAv
    step5 <- paste0(
      'Effective event probability after ', input$dr, '% total dropout:<br>',
      '<span style="font-family:monospace;">p\u2091 = p(event) \u00d7 (1 \u2212 ',
      sprintf("%.2f", input$dr/100), ') = ',
      sprintf("%.4f", pBase), ' \u00d7 ', sprintf("%.2f", 1-input$dr/100),
      ' = <strong>', sprintf("%.4f", r$peff), '</strong></span>'
    )
    
    step6 <- if (rv$arms == 1) paste0(
      'Sample size: N = \u2308D / p\u2091\u2309 = \u2308', r$D, ' / ', sprintf("%.4f", r$peff),
      '\u2309 = <strong style="font-family:monospace; color:', GN, ';">', r$N, '</strong>'
    ) else paste0(
      'Total N = \u2308D / p\u2091\u2309 = \u2308', r$Nt, ' / ', sprintf("%.4f", r$peff),
      '\u2309 = <strong style="font-family:monospace; color:', GN, ';">', r$Nt, '</strong><br>',
      'Treatment: n<sub>t</sub>  = \u2308N \u00b7 r\u2309 = <strong style="font-family:monospace; color:', TL, ';">', r$nt,
      '</strong> &nbsp;\u00b7&nbsp; Control: n<sub>c</sub> = N \u2212 n<sub>t</sub>  = <strong style="font-family:monospace; color:', OR, ';">', r$nc, '</strong>'
    )
    
    step7 <- paste0(
      'Verification: non-centrality \u03bc\u2081 = \u221aD \u00b7 |ln(HR)| = ', sprintf("%.4f", r$mu), '<br>',
      'Achieved power = \u03a6(\u03bc\u2081 \u2212 z\u2081\u208b\u03b1)',
      if(rv$sd==2) ' + \u03a6(\u2212\u03bc\u2081 \u2212 z\u2081\u208b\u03b1)' else '',
      ' = <span style="color:', GN, ';">', sprintf("%.2f%%", r$pw*100), '</span>',
      ' <span style="color:', DM, '; font-size:11px;">(target ',
      sprintf("%.0f%%", input$pw*100), '; small gap from \u2308D\u2309 ceiling rounding)</span>'
    )
    
    mk_li <- function(x) paste0('<li style="line-height:2.3;">', x, '</li>')
    HTML(paste0(
      '<div class="panel-box">',
      '<div class="sec-hdr">Step-by-step Calculation</div>',
      '<ol style="padding-left:18px; font-size:12px; color:', GY, ';">',
      mk_li(step1), mk_li(step2), mk_li(step3),
      mk_li(step4), mk_li(step5), mk_li(step6), mk_li(step7),
      '</ol></div>'
    ))
  })
  
  # ── Visualizations ────────────────────────────────────────────────────────
  output$svg_surv <- renderUI({
    req(is.finite(l0()), is.finite(l1()), input$L > input$A)
    HTML(make_surv_svg(l0(), l1(), input$A, input$L, rv$im, input$ts, rv$arms, U(), Ul()))
  })
  output$cap_surv <- renderText({
    if (rv$arms == 1)
      paste0("S(t) = exp(\u2212\u03bbt). H\u2080 benchmark (gray), H\u2081 target (orange). ",
             "Dashed at S = 0.5 with median drop-lines. ",
             "Pink = accrual A = ", input$A, " ", U(), ". Dashed vertical at L = ", input$L, " ", U(), ".")
    else
      paste0("S(t) = exp(\u2212\u03bbt). Control (gray), Treatment (orange). Dashed at S = 0.5. ",
             "Pink = accrual A = ", input$A, " ", U(), ". L = ", input$L, " ", U(), ".")
  })
  
  output$svg_z <- renderUI({
    if (!v()) return(NULL)
    HTML(make_z_svg(res()$mu, res()$z1a, rv$sd))
  })
  output$cap_z <- renderText({
    if (!v()) return("Enter valid inputs to display this plot.")
    r <- res()
    paste0("H\u2080: Z ~ N(0,1) gray. H\u2081: Z ~ N(\u03bc\u2081,1) blue, \u03bc\u2081 = ",
           sprintf("%.4f", r$mu), ". Gray = \u03b1 (false positive). Blue = \u03b2 (missed detection). ",
           "z(1\u2212\u03b1) = ", sprintf("%.4f", r$z1a), ". Power = ", sprintf("%.1f%%", r$pw*100), ".")
  })
  
  output$svg_evt <- renderUI({
    if (!v()) return(NULL)
    r <- res()
    N_val <- if (rv$arms == 1) r$N else r$Nt
    HTML(make_evt_svg(l1(), input$A, input$L, input$dr, N_val, l0(), rv$arms, input$k, U(), Ul()))
  })
  output$cap_evt <- renderText({
    if (!v()) return("Enter valid inputs to display this plot.")
    r <- res()
    N_val <- if (rv$arms == 1) r$N else r$Nt
    paste0("Projected cumulative events under H\u2081 with N = ", N_val, " enrolled, ",
           input$dr, "% total dropout, \u03bb\u2081 = ", sprintf("%.4f", l1()), " ", lU(),
           ", uniform accrual over A = ", input$A, " ", Ul(), ". ",
           "Pink dash = end of accrual. Green dash = total expected events at L = ", input$L, " ", U(), ".")
  })
  
  output$svg_heat <- renderUI({
    req(is.finite(l0()), is.finite(l1()), input$L > input$A)
    hr <- HR()
    cHR_val <- if (is.finite(hr)) hr else 0.7
    HTML(make_heat_svg(l0(), input$A, input$L, input$dr, input$al, rv$sd,
                       cHR_val, cN(), rv$arms, input$k))
  })
  output$cap_heat <- renderText({
    paste0("Power across sample size (N) and hazard ratio (HR = \u03bb\u2081/\u03bb\u2080) combinations, ",
           "holding \u03b1 = ", input$al, " (", if(rv$sd==1) "one" else "two", "-sided), ",
           "dropout = ", input$dr, "%, A = ", input$A, " ", U(), ", L = ", input$L, " ", U(), " fixed. ",
           "Orange = low power, blue = high power. White border = current design.")
  })
  
  # ── Sensitivity table ─────────────────────────────────────────────────────
  sens_rows <- reactive({
    req(is.finite(l0()), l0() > 0, is.finite(HR()))
    offsets <- -5:5
    lapply(offsets, function(off) {
      hv <- max(0.10, min(0.99, HR() + off * 0.05))
      ll  <- l0() * hv
      rowM1 <- ml_fn(ll)
      rowSt <- if (is.finite(ll) && input$ts > 0) Sv_fn(ll, input$ts) else NaN
      r <- if (rv$arms == 1)
        calcOne(l0(), ll, input$al, input$pw, rv$sd, input$A, input$L, input$dr)
      else
        calcTwo(l0(), ll, input$al, input$pw, rv$sd, rv$obj, input$dlt, input$k,
                input$A, input$L, input$dr)
      list(
        hv=round(hv,4), ll=ll, rowM1=rowM1, rowSt=rowSt,
        isActual=(off == 0),
        d  =if(!is.null(r)) r$D  else NaN,
        pe =if(!is.null(r)) (if(rv$arms==1) r$pe  else r$pAv)  else NaN,
        pef=if(!is.null(r)) r$peff else NaN,
        n  =if(!is.null(r)) (if(rv$arms==1) r$N   else r$Nt)   else NaN
      )
    })
  })
  
  output$sens_note <- renderUI({
    HTML(paste0(
      'Rows sweep \u00b15 steps of 0.05 HR around your current design (',
      switch(rv$im,
             "median" = paste0('m\u2080 = ', if(is.finite(m0())) sprintf("%.2f", m0()) else "\u2014", ' ', U(), ', m\u2081 varies'),
             "hazard" = paste0('\u03bb\u2080 = ', if(is.finite(l0())) sprintf("%.4f", l0()) else "\u2014", ' ', lU(), ', \u03bb\u2081 varies'),
             "St"     = paste0('S\u2080(t*=', input$ts, ') = ', sprintf("%.2f", input$s0), ', S\u2081 varies')
      ),
      '). \u03b1 = ', input$al, ' (', if(rv$sd==1) "one" else "two", '-sided), ',
      'power = ', sprintf("%.0f%%", input$pw*100), ', A = ', input$A, ' ', U(), ', L = ', input$L, ' ', U(), ', ',
      'dropout = ', input$dr, '%',
      if(rv$arms==2) paste0(', k = ', input$k) else '', '. ',
      '<strong style="color:', BL, ';">Middle row (blue) = your exact inputs.</strong>'
    ))
  })
  
  output$sens_table <- renderUI({
    rows <- sens_rows()
    if (length(rows) == 0) return(NULL)
    show_lam <- rv$im != "hazard"
    th <- function(txt, col=MU)
      paste0('<th style="padding:5px 9px; color:', col, '; text-align:left; ',
             'border-bottom:1px solid ', BDR, '; font-weight:500; font-size:11px; ',
             'white-space:nowrap;">', txt, '</th>')
    
    hdr <- paste0(
      '<tr>',
      th(switch(rv$im, "median"=paste0("m\u2080 (", U(), ") <span style='font-weight:400;'>H\u2080</span>"),
                "hazard"=paste0("\u03bb\u2080 (", lU(), ") <span style='font-weight:400;'>H\u2080</span>"),
                "St"    =paste0("S\u2080(t*=", input$ts, ") <span style='font-weight:400;'>H\u2080</span>")),
         GY),
      th(switch(rv$im, "median"=paste0("m\u2081 (", U(), ") <span style='font-weight:400;'>H\u2081</span>"),
                "hazard"=paste0("\u03bb\u2081 (", lU(), ") <span style='font-weight:400;'>H\u2081</span>"),
                "St"    =paste0("S\u2081(t*=", input$ts, ") <span style='font-weight:400;'>H\u2081</span>")),
         OR),
      th("HR = \u03bb\u2081/\u03bb\u2080"), th("ln(HR)"),
      if(show_lam) th(paste0("\u03bb\u2081 (", lU(), ")")) else "",
      th("Events (D)"),
      th(if(rv$arms==1) "P(event | H\u2081)" else "P\u0305(event)"),
      th("p\u2091 (dropout-adj.)"),
      th(if(rv$arms==1) "N" else "Total N", GN),
      '</tr>'
    )
    
    tds <- sapply(rows, function(row) {
      cur <- row$isActual
      h0Val <- switch(rv$im,
                      "median" = if(is.finite(m0())) sprintf("%.1f", m0()) else "\u2014",
                      "hazard" = if(is.finite(l0())) sprintf("%.4f", l0()) else "\u2014",
                      "St"     = sprintf("%.3f", input$s0)
      )
      h1Val <- switch(rv$im,
                      "median" = if(is.finite(row$rowM1)) sprintf("%.1f", row$rowM1) else "\u2014",
                      "hazard" = if(is.finite(row$ll))    sprintf("%.4f", row$ll)    else "\u2014",
                      "St"     = if(is.finite(row$rowSt)) sprintf("%.3f", row$rowSt) else "\u2014"
      )
      bg_row <- if(cur) paste0("background:#1a3a6a55; border-left:3px solid ", BL, ";") else
        "border-left:3px solid transparent;"
      td <- function(val, col=GY, fw=400, bgc="transparent")
        paste0('<td style="padding:5px 9px; font-family:monospace; color:', col,
               '; font-weight:', fw, '; background:', bgc, ';">', val, '</td>')
      paste0('<tr style="', bg_row, '">',
             td(h0Val, GY),
             td(h1Val, if(cur) OR else MU, if(cur) 700 else 400, if(cur) paste0(OR,"22") else "transparent"),
             td(sprintf("%.4f", row$hv),   if(cur) BL else GY, if(cur) 700 else 400),
             td(sprintf("%.3f", log(row$hv)), MU),
             if(show_lam) td(if(is.finite(row$ll)) sprintf("%.4f", row$ll) else "\u2014", MU) else "",
             td(if(is.finite(row$d)) as.integer(row$d) else "\u2014", if(cur) BL else TX, if(cur) 700 else 400),
             td(if(is.finite(row$pe)) sprintf("%.1f%%", row$pe*100) else "\u2014", MU),
             td(if(is.finite(row$pef)) sprintf("%.1f%%", row$pef*100) else "\u2014", MU),
             td(if(is.finite(row$n)) as.integer(row$n) else "\u2014", if(cur) GN else TX, if(cur) 700 else 400),
             '</tr>'
      )
    })
    
    HTML(paste0(
      '<table class="sens-table"><thead>', hdr, '</thead><tbody>',
      paste(tds, collapse=""), '</tbody></table>'
    ))
  })
  
  # ── Budget back-calculation ───────────────────────────────────────────────
  bc <- reactive({
    n <- input$budget_n
    if (is.null(n) || is.na(n) || n < 1) return(NULL)
    if (!is.finite(l0()) || !is.finite(l1())) return(NULL)
    backCalc(round(n), l0(), l1(), input$al, input$pw, rv$sd, rv$obj,
             input$dlt, input$k, input$A, input$L, input$dr, rv$arms)
  })
  
  output$budget_result <- renderUI({
    b <- bc()
    if (is.null(b)) return(NULL)
    n <- round(input$budget_n)
    
    if (!b$feasible) {
      return(HTML(paste0(
        '<div style="background:#2d0a0a; border:1px solid ', RD, ';',
        'border-radius:8px; padding:11px 15px; color:', RD, '; font-size:13px;">',
        '<strong>Not feasible:</strong> ', b$reason, '</div>'
      )))
    }
    
    pw_col <- if (b$power >= 0.80) GN else if (b$power >= 0.60) OR else RD
    card <- function(lbl, val, sub, col)
      paste0('<div style="background:', BG, '; border:1px solid ', col, '33;',
             'border-radius:8px; padding:10px 14px; flex:1 1 100px;">',
             '<div style="color:', MU, '; font-size:10px; letter-spacing:1px;',
             'text-transform:uppercase; margin-bottom:4px;">', lbl, '</div>',
             '<div style="color:', col, '; font-size:20px; font-weight:700;',
             'font-family:monospace;">', val, '</div>',
             '<div style="color:', MU, '; font-size:11px; margin-top:3px;">', sub, '</div>',
             '</div>')
    
    alpha_card <- if (!b$impossible && !is.null(b$alpha1s)) {
      card(paste0("Alpha needed for ", sprintf("%.0f%%", input$pw*100), "power"), #issue: \u03b1 not rendering correctly
           sprintf("%.4f", b$alpha1s),
           paste0("one-sided",
                  if(!is.null(b$alpha2s)) paste0(" (", sprintf("%.4f", b$alpha2s), " two-sided)") else ""),
           OR)
    } else ""
    
    cards <- paste0(
      '<div style="display:flex; gap:10px; flex-wrap:wrap; margin-bottom:14px;">',
      card("Expected Events (D)", b$D,
           paste0("from N = ", n, " \u00d7 p\u2091 = ", sprintf("%.1f%%", b$pe*100)),
           BL),
      card("Achieved Power",
           paste0(sprintf("%.1f%%", b$power*100)),
           paste0("at \u03b1 = ", input$al, " (", if(rv$sd==1) "one" else "two", "-sided)"),
           pw_col),
      alpha_card,
      '</div>'
    )
    
    interp <- if (b$impossible) {
      paste0('<strong style="color:', RD, ';">Target power not achievable.</strong> ',
             'With N = ', n, ', the expected ', b$D, ' event', if(b$D!=1) "s" else "",
             ' yield a non-centrality of \u03bc\u2081 = ', sprintf("%.3f", b$mu1),
             ', which is less than z\u2081\u208b\u03b2 = ', sprintf("%.3f", qnorm(input$pw)),
             ' required for ', sprintf("%.0f%%", input$pw*100), 'power. ',
             'To reach ', sprintf("%.0f%%", input$pw*100), 'power you would need to ',
             'increase N, extend follow-up, or accept a larger HR difference.')
    } else {
      ai <- if (!is.null(b$alpha1s)) {
        commentary <- if (b$alpha1s > 0.10) " This \u03b1 is unusually large and may not be scientifically justifiable."
        else if (b$alpha1s > 0.05) " This is a modest relaxation of the standard \u03b1 = 0.05."
        else " This is within conventional bounds."
        paste0(' To achieve your target of ', sprintf("%.0f%%", input$pw*100), 'power, ',
               'the one-sided \u03b1 would need to be relaxed to ',
               '<strong style="color:', OR, ';">', sprintf("%.4f", b$alpha1s), '</strong>',
               if(!is.null(b$alpha2s)) paste0(' (', sprintf("%.4f", b$alpha2s), ' two-sided)') else "",
               '.', commentary)
      } else ""
      paste0('With a budget of <strong style="color:', BL, ';">N = ', n, '</strong>, ',
             'approximately <strong style="color:', BL, ';">', b$D, ' events</strong> ',
             'are expected (p\u2091 = ', sprintf("%.1f%%", b$pe*100), '%). ',
             'At your original \u03b1 = ', input$al, ' (',
             if(rv$sd==1) "one" else "two", '-sided), this yields a power of ',
             '<strong style="color:', pw_col, ';">', sprintf("%.1f%%", b$power*100), '</strong>.', ai)
    }
    
    HTML(paste0(
      cards,
      '<div style="background:#0d1f38; border:1px solid ', B2, ';',
      'border-radius:8px; padding:11px 15px; font-size:13px; color:', TX, '; line-height:1.9;">',
      interp, '</div>'
    ))
  })
  
  # ── References ────────────────────────────────────────────────────────────
  output$references <- renderUI({
    refs <- list(
      list(cite="Schoenfeld DA (1981).",
           title="The asymptotic properties of nonparametric tests for comparing survival distributions.",
           journal="Biometrika 68(1):316\u2013319.",
           url="https://doi.org/10.1093/biomet/68.1.316"),
      list(cite="Schoenfeld DA (1983).",
           title="Sample-size formula for the proportional-hazards regression model.",
           journal="Biometrics 39(2):499\u2013503.",
           url="https://doi.org/10.2307/2531021"),
      list(cite="Lachin JM, Foulkes MA (1986).",
           title="Evaluation of sample size and power for analyses of survival with allowance for nonuniform patient entry, losses to follow-up, noncompliance, and stratification.",
           journal="Biometrics 42(3):507\u2013519.",
           url="https://doi.org/10.2307/2531201"),
      list(cite="Harrington DP, Fleming TR (1982).",
           title="A class of rank test procedures for censored survival data.",
           journal="Biometrika 69(3):553\u2013566.",
           url="https://doi.org/10.1093/biomet/69.3.553"),
      list(cite="Jung S-H, Chow S-C, Chi EM (2005).",
           title="A note on sample size calculation based on propensity analysis in non-randomized trials.",
           journal="Biometrics 61(3):844\u2013848.",
           url="https://doi.org/10.1111/j.1541-0420.2005.00345.x"),
      list(cite="Zhou Y, Yuan Y, Lee JJ, Pan H.",
           title="Nsurvival: Sample size calculation for time-to-event endpoints.",
           journal="MD Anderson Cancer Center Shiny App.",
           url="https://biostatistics.mdanderson.org/shinyapps/Nsurvival/")
    )
    items <- sapply(seq_along(refs), function(i) {
      r <- refs[[i]]
      sep <- if (i < length(refs)) paste0("border-bottom:1px solid ", BDR, ";") else ""
      paste0(
        '<div style="display:flex; gap:10px; align-items:flex-start; padding:7px 0; ', sep, '">',
        '<span style="color:', B2, '; font-size:11px; flex-shrink:0; font-family:monospace; padding-top:1px;">[', i, ']</span>',
        '<span style="line-height:1.7;">',
        '<span style="color:', GY, ';">', r$cite, ' </span>',
        '<span style="color:', TX, ';">', r$title, ' </span>',
        '<span style="color:', MU, '; font-style:italic;">', r$journal, ' </span>',
        '<a href="', r$url, '" target="_blank" rel="noopener noreferrer" ',
        'style="color:', BL, '; text-decoration:none; font-size:10px;">', r$url, '</a>',
        '</span></div>'
      )
    })
    HTML(paste0('<div style="font-size:11px; color:', GY, '; line-height:1.5;">',
                paste(items, collapse=""), '</div>'))
  })
}

shinyApp(ui, server)
