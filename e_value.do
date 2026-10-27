* Author: Yang Liu
* Date: 2026-10-27
* Stata version: 18 (compatible with 14 or later)
*
* Prepared input ----
* do e_value.do "data/inputs/cox_model3.csv" "data/outputs/e_values.csv"
* Columns: exposure outcome contrast hr lcl ucl participants cases.
* Fully adjusted, full-precision HRs: five outcomes x three contrasts.
* Uses the rare-outcome HR-as-RR approximation.

version 14
* E-values ----
args input_csv output_csv
if `"`input_csv'"' == "" {
    display as error "Supply input CSV and output CSV as arguments."
    exit 198
}
if `"`output_csv'"' == "" {
    display as error "Supply an output CSV path; its parent directory must exist."
    exit 198
}
preserve
import delimited using `"`input_csv'"', varnames(1) clear case(preserve)
confirm variable exposure outcome contrast hr lcl ucl participants cases
assert _N == 15
assert exposure == "Benzene"
assert inlist(outcome, "All_CB", "AB", "LB", "RB", "Other_CB")
assert inlist(contrast, "T2 vs T1", "T3 vs T1", "Per SD increase")
isid exposure outcome contrast
bysort exposure outcome: assert _N == 3
assert !missing(hr,lcl,ucl,participants,cases)
assert hr > 0 & lcl > 0 & lcl <= hr & hr <= ucl
assert participants > 0 & cases >= 0 & cases <= participants
gen double observed_event_fraction = cases / participants
assert observed_event_fraction < 0.15
gen double rr_point = cond(hr >= 1, hr, 1/hr)
gen double rr_ci = 1
replace rr_ci = lcl if hr > 1 & lcl > 1
replace rr_ci = 1/ucl if hr < 1 & ucl < 1
gen double evalue_point = rr_point + sqrt(rr_point*(rr_point-1))
gen double evalue_ci = rr_ci + sqrt(rr_ci*(rr_ci-1))
assert evalue_point >= 1 & evalue_ci >= 1 & evalue_ci <= evalue_point
format hr lcl ucl evalue_point evalue_ci %24.17g
export delimited exposure outcome contrast hr lcl ucl participants cases ///
    observed_event_fraction evalue_point evalue_ci using `"`output_csv'"', replace
restore
