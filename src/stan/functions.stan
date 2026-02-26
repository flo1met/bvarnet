functions {
    real set_prior(vector x, int fam, real loc, real scale, real df) {
        if (fam == 1) return normal_lpdf(x | loc, scale);
        if (fam == 2) return student_t_lpdf(x | df, loc, scale);
        if (fam == 3) return cauchy_lpdf(x | loc, scale);
        reject("Unknown prior family code:", fam);
    }

    real set_half_prior(vector x, int fam, real loc, real scale, real df) {
        if (fam == 1) return normal_lpdf(x | loc, scale) - rows(x) * normal_lccdf(0 | loc, scale);
        if (fam == 2) return student_t_lpdf(x | df, loc, scale) - rows(x) * student_t_lccdf(0 | df, loc, scale);
        if (fam == 3) return cauchy_lpdf(x | loc, scale) - rows(x) * cauchy_lccdf(0 | loc, scale);
        reject("Unknown prior family code:", fam);
   }
}