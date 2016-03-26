A = pow(kT,2)*(pow(kT,2) + pow(pT,2) + pow(qT,2) - 
      2*kT*pT*cos(phi) + 2*pT*qT*cos(phiq) - 
      2*kT*qT*cos(phi + phiq))*
    (1/(4.*pow(pT,2)*pow(qT,2)) - 
      pow(E,dyQMRK)/
       (4.*pT*qT*(-pow(kT,2) - pow(pT,2) - 
           pow(E,dyQMRK)*pT*qT + 2*kT*pT*cos(phi))) - 
      1/(4.*pow(E,dyQMRK)*pT*qT*
         (-pow(kT,2) - (pT*qT)/pow(E,dyQMRK) - pow(qT,2) + 
           2*kT*qT*cos(phi + phiq))) - 
      1/((-pow(kT,2) - pow(pT,2) - pow(E,dyQMRK)*pT*qT + 
           2*kT*pT*cos(phi))*
         (-pow(kT,2) - (pT*qT)/pow(E,dyQMRK) - pow(qT,2) + 
           2*kT*qT*cos(phi + phiq))) + 
      (pow(kT,2)*(pow(kT,2) + pow(pT,2) + pow(qT,2) - 
           2*kT*pT*cos(phi) + 2*pT*qT*cos(phiq) - 
           2*kT*qT*cos(phi + phiq)))/
       (4.*pow(pT,2)*pow(qT,2)*
         (-pow(kT,2) - pow(pT,2) - pow(E,dyQMRK)*pT*qT + 
           2*kT*pT*cos(phi))*
         (-pow(kT,2) - (pT*qT)/pow(E,dyQMRK) - pow(qT,2) + 
           2*kT*qT*cos(phi + phiq))) + 
      (-(pow(kT,2)*((1 + qT/(pow(E,dyQMRK)*pT))/
                (-pow(kT,2) - pow(pT,2) - pow(E,dyQMRK)*pT*qT + 
                  2*kT*pT*cos(phi)) + 
               (1 + (pow(E,dyQMRK)*pT)/qT)/
                (-pow(kT,2) - (pT*qT)/pow(E,dyQMRK) - 
                  pow(qT,2) + 2*kT*qT*cos(phi + phiq))))/
          (8.*(-(pT*qT*cos(phiq)) + pT*qT*cosh(dyQMRK))) - 
         ((pow(kT,2) + pow(pT,2) + pow(qT,2) - 
              2*kT*pT*cos(phi) + 2*pT*qT*cos(phiq) - 
              2*kT*qT*cos(phi + phiq))*
            ((1 + pT/(pow(E,dyQMRK)*qT))/
               (-pow(kT,2) - pow(pT,2) - pow(E,dyQMRK)*pT*qT + 
                 2*kT*pT*cos(phi)) + 
              (1 + (pow(E,dyQMRK)*qT)/pT)/
               (-pow(kT,2) - (pT*qT)/pow(E,dyQMRK) - 
                 pow(qT,2) + 2*kT*qT*cos(phi + phiq))))/
          (8.*(-(pT*qT*cos(phiq)) + pT*qT*cosh(dyQMRK))) + 
         (cosh(dyQMRK)*(1 + (pow(pT,2) + pow(qT,2) + 
                 2*pT*qT*cosh(dyQMRK))/
               (2.*(-(pT*qT*cos(phiq)) + pT*qT*cosh(dyQMRK)))))/
          (2.*pT*qT) - (1 + 
            pT*qT*(1/
                (-pow(kT,2) - pow(pT,2) - pow(E,dyQMRK)*pT*qT + 
                  2*kT*pT*cos(phi)) - 
               1/
                (-pow(kT,2) - (pT*qT)/pow(E,dyQMRK) - 
                  pow(qT,2) + 2*kT*qT*cos(phi + phiq)))*
             sinh(dyQMRK))/(-(pT*qT*cos(phiq)) + pT*qT*cosh(dyQMRK)))/
       (pow(pT,2) + pow(qT,2) + 2*pT*qT*cosh(dyQMRK))) + 
   (pow(-(pow(pT,2)*pow(qT,2)) + 
         (pow(kT,2) + pow(pT,2) - 2*kT*pT*cos(phi))*
          (pow(kT,2) + pow(qT,2) - 2*kT*qT*cos(phi + phiq)),
        2)/(pow(-pow(kT,2) - pow(pT,2) - 
           pow(E,dyQMRK)*pT*qT + 2*kT*pT*cos(phi),2)*
         pow(-pow(kT,2) - (pT*qT)/pow(E,dyQMRK) - pow(qT,2) + 
           2*kT*qT*cos(phi + phiq),2)) - 
      (((pow(kT,2) - (pT*qT)/pow(E,dyQMRK) + pow(qT,2) - 
              2*kT*qT*cos(phi + phiq))/
            (pow(kT,2) + (pT*qT)/pow(E,dyQMRK) + pow(qT,2) - 
              2*kT*qT*cos(phi + phiq)) - 
           (-pow(pT,2) + pow(qT,2) + 2*kT*pT*cos(phi) - 
              2*kT*qT*cos(phi + phiq) - 
              ((pow(pT,2) - pow(qT,2))*
                 (-pow(pT,2) - pow(qT,2) + 2*kT*pT*cos(phi) - 
                   2*pT*qT*cos(phiq) + 2*kT*qT*cos(phi + phiq)))
                /(pow(pT,2) + pow(qT,2) + 2*pT*qT*cosh(dyQMRK)) + 
              2*pT*qT*(1 - 
                 (2*pow(kT,2) + pow(pT,2) + pow(qT,2) - 
                    2*kT*pT*cos(phi) + 2*pT*qT*cos(phiq) - 
                    2*kT*qT*cos(phi + phiq))/
                  (pow(pT,2) + pow(qT,2) + 2*pT*qT*cosh(dyQMRK)))*
               sinh(dyQMRK))/
            (2.*(-(pT*qT*cos(phiq)) + pT*qT*cosh(dyQMRK))))*
         ((pow(kT,2) + pow(pT,2) - pow(E,dyQMRK)*pT*qT - 
              2*kT*pT*cos(phi))/
            (pow(kT,2) + pow(pT,2) + pow(E,dyQMRK)*pT*qT - 
              2*kT*pT*cos(phi)) + 
           (-pow(pT,2) + pow(qT,2) + 2*kT*pT*cos(phi) - 
              2*kT*qT*cos(phi + phiq) - 
              ((pow(pT,2) - pow(qT,2))*
                 (-pow(pT,2) - pow(qT,2) + 2*kT*pT*cos(phi) - 
                   2*pT*qT*cos(phiq) + 2*kT*qT*cos(phi + phiq)))
                /(pow(pT,2) + pow(qT,2) + 2*pT*qT*cosh(dyQMRK)) + 
              2*pT*qT*(1 - 
                 (2*pow(kT,2) + pow(pT,2) + pow(qT,2) - 
                    2*kT*pT*cos(phi) + 2*pT*qT*cos(phiq) - 
                    2*kT*qT*cos(phi + phiq))/
                  (pow(pT,2) + pow(qT,2) + 2*pT*qT*cosh(dyQMRK)))*
               sinh(dyQMRK))/
            (2.*(-(pT*qT*cos(phiq)) + pT*qT*cosh(dyQMRK)))))/4.)/2.;
