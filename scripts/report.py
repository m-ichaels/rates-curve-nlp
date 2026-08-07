#!/usr/bin/env python3
"""report.pdf from results/summary.md and results/figures/*.png (fpdf2).   python scripts/report.py"""
import os
import re
import sys

from fpdf import FPDF

R = sys.argv[1] if len(sys.argv) > 1 else "results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "report.pdf"

INTRO = """Question. A Treasury RFQ desk needs a curve to price, a curve to hedge and a curve to quote off. Does machine learning improve on the classical answers - a Hagan-West bootstrap, key-rate DV01 hedges, closed-form quotes (Bergault et al., Barzykin-Ciceri, Cartea-Wang) - when each is tested out of sample on 2015-2026 Treasury data and, for quoting, in an RFQ simulator with informed flow and scheduled events?

Method. Header-only C++ behind one CLI. Curve: monotone-convex bootstrap of the Treasury CMT curve with closed-form sector integrals, key-rate DV01s, a meeting-date step front end, and a kernel-ridge discount curve (Filipovic-Pelger-Ye) with the smoothing chosen by a rule fixed in advance. Daily builds 2015-2026 scored by leave-one-out pricing error, forward roughness, DV01 stability and a hedge test (an off-curve bond hedged with the curve's own key-rate DV01s, realised next-day P&L). Hedge ratios learned from co-movements (ridge, PCA) tested walk-forward against the model ratios. Risk: PCA factors; RFQ simulator with tiered clients, informed pre-event flow, event jumps, band hedging, exact P&L decomposition. Quoters: the closed forms and an oracle closed form that knows the client model. Events: FOMC/ECB statement hawkishness (lexicon; FOMC-RoBERTa dated checkpoint), USMPD reactions, NFP/auction calendar, CFTC positioning, funding stress; text signal fitted on 2015-2023 and replayed on 2024-2026 with permutation and look-ahead placebos.

Caveats. SOFR OIS and CME futures are not free: the curve is the Treasury curve and the step front end is fitted from bills. Client flow is simulated. Few meetings per year: text CIs are wide by construction."""

FIGS = [("curve.png", "Curve construction on one day: bootstrap vs kernel-ridge forwards and zeros, leave-one-out pricing error (interior tenors), meeting-date step front end."),
        ("curvehist.png", "Daily builds 2015-2026: 10y par DV01 under both methods, leave-one-out error by year, hedge test with each curve's key-rate DV01s, smoothness vs error."),
        ("hedge.png", "Learned vs model hedge ratios, walk-forward: ridge on the trailing window vs the curves' key-rate DV01s, PCA-neutral and 50/50."),
        ("factors.png", "PCA loadings by window: level, slope, curvature and their daily standard deviations."),
        ("frontier.png", "Closed-form quoters on synthetic days: frontier over risk aversion, P&L decomposition, mark-outs by tier around events."),
        ("policy.png", "Quote structure of the closed forms: mean half-spread by tier and by risk-reducing / risk-adding side, normal RFQs and the pre-event window."),
        ("nlp.png", "NLP arm: out-of-sample correlation of the text nowcast with the event-window move by tenor (ticks = permutation placebo), 2y nowcast scatter, dated-checkpoint split, size-of-move model."),
        ("text_signal.png", "Text signal: change in hawkishness vs 10y reaction (train/test), hawkishness through time, look-ahead placebo."),
        ("events.png", "Real-day replay 2024-2026 with the event features: P&L, P&L per unit variance with signal / permuted signal / no signal, mark-outs around events.")]


class PDF(FPDF):
    def header(self):
        self.set_font("Helvetica", "B", 9); self.set_text_color(120); self.cell(0, 6, "rates-curve-nlp - machine-learned yield curves for pricing, hedging and quoting", align="R"); self.ln(8); self.set_text_color(0)

    def footer(self):
        self.set_y(-12); self.set_font("Helvetica", "", 8); self.set_text_color(120); self.cell(0, 6, f"{self.page_no()}", align="C")


def clean(s):
    return (s.replace("–", "-").replace("—", "-").replace("−", "-").replace("σ", "sigma").replace("β", "beta").replace("η", "eta").replace("Δ", "d").replace("×", "x").replace("≥", ">=").replace("…", "...")
             .replace("²", "^2").replace("±", "+/-").replace("**", "").replace("`", "").replace("√", "sqrt ").replace("γ", "gamma").replace("ψ", "psi").replace("μ", "mu").replace("ρ", "rho").replace("λ", "lambda").replace("→", "->").replace("≈", "~"))


def md_table(pdf, rows):
    cols = [c.strip() for c in rows[0].strip("|").split("|")]
    data = [[clean(c.strip()) for c in r.strip("|").split("|")] for r in rows[2:]]
    pdf.set_font("Helvetica", "", 6.5)
    n = len(cols); w = (pdf.w - 20) / n
    widths = [w] * n
    pdf.set_font("Helvetica", "B", 6.5)
    for c, wd in zip(cols, widths): pdf.cell(wd, 5, clean(c)[:40], border=1)
    pdf.ln(5); pdf.set_font("Helvetica", "", 6.5)
    for r in data:
        if pdf.get_y() > pdf.h - 20: pdf.add_page()
        for c, wd in zip(r, widths): pdf.cell(wd, 4.5, c[:40], border=1)
        pdf.ln(4.5)
    pdf.ln(2)


def main():
    pdf = PDF(); pdf.set_auto_page_break(auto=True, margin=15); pdf.add_page()
    pdf.set_font("Helvetica", "B", 16); pdf.cell(0, 10, "Machine-Learned Yield Curves for Pricing, Hedging and Quoting", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 9)
    for para in INTRO.split("\n\n"):
        pdf.multi_cell(0, 4.5, clean(para)); pdf.ln(2)
    for fn, cap in FIGS:
        p = os.path.join(R, "figures", fn)
        if not os.path.exists(p): continue
        if pdf.get_y() > pdf.h - 90: pdf.add_page()
        pdf.image(p, w=pdf.w - 20); pdf.set_font("Helvetica", "I", 8); pdf.multi_cell(0, 4, clean(cap)); pdf.ln(3); pdf.set_font("Helvetica", "", 9)
    # summary tables
    sm = os.path.join(R, "summary.md")
    if os.path.exists(sm):
        pdf.add_page()
        lines = open(sm, encoding="utf-8").read().splitlines()
        i = 0
        while i < len(lines):
            l = lines[i]
            if l.startswith("## "):
                pdf.set_font("Helvetica", "B", 11); pdf.ln(2); pdf.cell(0, 7, clean(l[3:]), new_x="LMARGIN", new_y="NEXT"); pdf.set_font("Helvetica", "", 9); i += 1
            elif l.startswith("### "):
                pdf.set_font("Helvetica", "B", 9); pdf.cell(0, 6, clean(l[4:]), new_x="LMARGIN", new_y="NEXT"); pdf.set_font("Helvetica", "", 9); i += 1
            elif l.startswith("|"):
                j = i
                while j < len(lines) and lines[j].startswith("|"): j += 1
                if j - i >= 2: md_table(pdf, lines[i:j])
                i = j
            elif l.startswith("```"):
                j = i + 1
                while j < len(lines) and not lines[j].startswith("```"): j += 1
                pdf.set_font("Courier", "", 7)
                for t in lines[i + 1:j]: pdf.set_x(pdf.l_margin); pdf.multi_cell(0, 3.5, clean(t))
                pdf.set_font("Helvetica", "", 9); i = j + 1
            elif l.strip():
                pdf.set_x(pdf.l_margin); pdf.multi_cell(0, 4.5, clean(l.strip())); i += 1
            else:
                i += 1
    pdf.output(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
