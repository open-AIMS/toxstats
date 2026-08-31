# the Appendix C summary is stable

    Code
      summary(toxcalc(fathead_c1, response = "weight", pmsd_bounds = "fathead_growth"))
    Output
      EPA WET hypothesis test
      Flowchart: EPA-821-R-02-013 Figure 2
      
         1  Is the response a proportion requiring transformation?
            no transformation needed
            (EPA-821-R-02-013 Appendix B, section 4.2)
         2  Are the pooled within-group residuals normally distributed?
            0.9507, p = 0.378 -> residuals consistent with normality
            (EPA-821-R-02-013 Appendix B, section 2.1)
         3  Are the variances homogeneous across concentrations?
            7.856, p = 0.097 -> variances not significantly different
            (EPA-821-R-02-013 Appendix B, section 3)
         4  Is replication equal across all concentrations?
            4, 4, 4, 4, 4 replicates, balanced
            (EPA-821-R-02-013 Figure 2)
         5  Which test was run?
            Dunnett's procedure
            (EPA-821-R-02-013 Appendix C)
      
        NOEC 128    LOEC 256
        MSD  0.1618    PMSD 23.9 per cent (EPA bounds 12 to 30: within)

# the Appendix E summary is stable

    Code
      summary(toxcalc(ceriodaphnia_e1, response = "young", exclude = 50))
    Output
      EPA WET hypothesis test
      Flowchart: EPA-821-R-02-013 Figure 2
      
      Excluded from the hypothesis test:
        50 -- excluded by the caller (section 9.5.2)
      
         1  Is the response a proportion requiring transformation?
            no transformation needed
            (EPA-821-R-02-013 Appendix B, section 4.2)
         2  Are the pooled within-group residuals normally distributed?
            0.9282, p = 0.00471 -> residuals not normally distributed
            (EPA-821-R-02-013 Appendix B, section 2.1)
         3  Are there at least four replicates at every concentration?
            10 -> smallest concentration has 10 replicates
            (EPA-821-R-02-013 section 9.4.5.2)
         4  Is replication equal across all concentrations?
            10, 10, 10, 10, 10 replicates, balanced
            (EPA-821-R-02-013 Figure 2)
         5  Which test was run?
            Steel's Many-One Rank Test (asymptotic p-values)
            (EPA-821-R-02-013 Appendix E)
      
        NOEC 3    LOEC 6
        MSD  4.92    PMSD 20 per cent

