#!/bin/sh

if [ -e curves_analytic/                   ]; then rm -r curves_analytic/                   ; fi
if [ -e curves_analytic_unfolded/          ]; then rm -r curves_analytic_unfolded/          ; fi
if [ -e curves_custom/                     ]; then rm -r curves_custom/                     ; fi
if [ -e curves_gmct/                       ]; then rm -r curves_gmct/                       ; fi
if [ -e page01.gnuplot                     ]; then rm -r page01.gnuplot                     ; fi
if [ -e page02.gnuplot                     ]; then rm -r page02.gnuplot                     ; fi
if [ -e page03.gnuplot                     ]; then rm -r page03.gnuplot                     ; fi
if [ -e page04.gnuplot                     ]; then rm -r page04.gnuplot                     ; fi
if [ -e prob_ph7_analytic.sed              ]; then rm -r prob_ph7_analytic.sed              ; fi
if [ -e prob_ph7_analytic_unfolded.sed     ]; then rm -r prob_ph7_analytic_unfolded.sed     ; fi
if [ -e prob_ph7_custom.sed                ]; then rm -r prob_ph7_custom.sed                ; fi
if [ -e prob_ph7_gmct.sed                  ]; then rm -r prob_ph7_gmct.sed                  ; fi
if [ -e substate_ph7_analytic.tex          ]; then rm -r substate_ph7_analytic.tex          ; fi
if [ -e substate_ph7_analytic_unfolded.tex ]; then rm -r substate_ph7_analytic_unfolded.tex ; fi
if [ -e substate_ph7_custom.tex            ]; then rm -r substate_ph7_custom.tex            ; fi
if [ -e substate_ph7_gmct.tex              ]; then rm -r substate_ph7_gmct.tex              ; fi
if [ -e rubredoxin.out                     ]; then rm -r rubredoxin.out                     ; fi
if [ -e plots.pdf                          ]; then rm -r plots.pdf                          ; fi
if [ -e mead/                              ]; then rm -r mead/                              ; fi
