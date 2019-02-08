Feature: formula is extracted properly from Basis class.

  Scenario:
    When Basis is created with the following data:
      | basis                               | cacheKey |
      | $BASIS{Si 0 0 0, Ge 0.25 0.25 0.25} | basis    |
    Then Basis stored under "basis" key returns "SiGe" as formula
    
