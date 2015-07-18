#!/bin/sh

if [ -e curves_mc/       ]; then rm -r curves_mc/       ; fi
if [ -e curves_analytic/ ]; then rm -r curves_analytic/ ; fi
if [ -e page01.gnuplot   ]; then rm -r page01.gnuplot   ; fi
if [ -e page02.gnuplot   ]; then rm -r page02.gnuplot   ; fi
if [ -e page03.gnuplot   ]; then rm -r page03.gnuplot   ; fi
if [ -e page04.gnuplot   ]; then rm -r page04.gnuplot   ; fi
if [ -e lysozyme.out     ]; then rm -r lysozyme.out     ; fi
if [ -e plots.pdf        ]; then rm -r plots.pdf        ; fi
if [ -e mead/            ]; then rm -r mead/            ; fi
