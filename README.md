This is an example to run the main function ConvFunc.


Assuming that FSIR analysis for a binary trait with a prevalence K=0.3 was performed and the following estimates on the observed 0-1 scale were obtained: hsq = 0.224 (standard error (s.e.) 0.0189) and csq = 0.0999 (s.e. 0.0374).

The code below takes approximately 10s to generate an vector of converted estimates on the liability scale.

`results = ConvFun(hsq_01=0.224,se_hsq_01=0.0189,csq_01=0.0999,se_csq_01=0.0374,K=0.3)`
`print(results)`
`       hsq_l     se_hsq_l        csq_l     se_csq_l         Crit 
3.400342e-01 2.806736e-02 1.764542e-01 6.205887e-02 4.079718e-09 `
