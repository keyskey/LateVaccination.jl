## Presumption: This is perfect immunity model

### State change
States are;
1. IM
2. S
3. I
4. R

- Pre-emptive vaccinators(PV) -> Immuned(Doesn't change state)  
- Free rider(S) -> Immuned (Went to hospital, Late vaccinator(LV))  
                -> Infected -> Recovered (Didn't go to hospital and infected, failed free rider)    
                -> Free rider (Keep S, successful free rider)

P(Free rider -> Immunued) = num_i/num_im  

### Strategy update
- IBRA when changing season  
- Decide whether go to Pre-emptive vaccinators or Late vaccinators
