
/* A solvent with the three Hansen parameters. */
struct Solvent {
    d_param: f64,
    p_param: f64,
    h_param: f64
}

impl Solvent {
    
    fn scale(&self, factor: f64) -> Solvent {
        Solvent {
            d_param: self.d_param * factor,
            p_param: self.p_param * factor,
            h_param: self.h_param * factor}
    }
    
    fn add(&self, other: Solvent) -> Solvent {
        Solvent {
            d_param: self.d_param + other.d_param,
            p_param: self.p_param + other.p_param,
            h_param: self.h_param + other.h_param}
    }
    
    fn r_value(&self, other: Solvent) -> f64 {
        let diffs = vec![
            self.d_param - other.d_param,
            self.p_param - other.p_param,
            self.h_param - other.h_param];
        diffs.iter()
            .map(|&param| param * param)
            .fold(0.0f64, |sum, param| sum + param)
    }
    
    fn combine(proportionated: &[(Solvent, f64)]) -> Solvent {
        let sum_proportion =
            proportionated.iter()
                .map(|&sp| { let (_,p) = sp; p })
                .fold(0.0f64, |sum, p| sum + p );
        
        proportionated.iter()
            .map( |&sp| { let (s,p) = sp; s.scale(p) } )
            .fold(Solvent::null_solvent(), |s1,s2| s1.add(s2))
            .scale(1.0f64/sum_proportion)
    }
    
    /* Because monoid. */
    fn null_solvent() -> Solvent {
        Solvent { d_param: 0.0, p_param: 0.0, h_param: 0.0 }
    }
    
    fn acetic_acid() -> Solvent {
        Solvent { d_param: 6.8, p_param: 6.0, h_param: 9.2 }
    }
    
    fn dioctyl_phtalate() -> Solvent {
        Solvent { d_param: 8.1, p_param: 3.4, h_param: 1.5 }
    }
    
    fn nitroamine() -> Solvent {
        Solvent { d_param: 8.6, p_param: 6.8, h_param: 0.0 }
    }
    
    /* Mystery solvent MEP */
    fn mep() -> Solvent {
        Solvent { d_param: 7.8, p_param: 4.4, h_param: 2.5 }
    }
    
}

fn main() {
    let prospective =
        Solvent::combine([
            (Solvent::acetic_acid(),      2.0),
            (Solvent::dioctyl_phtalate(), 5.0),
            (Solvent::nitroamine(),       3.0)
        ]);
    
    println!("R-value for 20% AA, 50% DP 30% NA: {}", 
            prospective.r_value(Solvent::mep()));
    
    let mut best_proportion = (0f64, 0f64, 0f64);
    let mut best_r_value = 1000f64;
    let target_solvent = Solvent::mep();
    
    let precision = 100i;
    
    for aa in std::iter::range_inclusive(0i, precision) {
        for dp in std::iter::range_inclusive(0i, precision) {
            if aa + dp <= precision {
                let na = precision - aa - dp;
                let prospective =
                    Solvent::combine([
                        (Solvent::acetic_acid(),      aa as f64),
                        (Solvent::dioctyl_phtalate(), dp as f64),
                        (Solvent::nitroamine(),       na as f64)
                    ]);
                
                let current_r = prospective.r_value(target_solvent);
                
                if current_r < best_r_value {
                    best_proportion = (aa as f64, dp as f64, na as f64);
                    best_r_value = current_r;
                }
            }
        }
    }
    
    let (aa, dp, na) = best_proportion;
    let sum = aa + dp + na;
    println!("Best match for MEP is:");
    println!("{}% Acetic Acid", aa/sum*100.0);
    println!("{}% Dioctyl Phtalate", dp/sum*100.0);
    println!("{}% Nitroamine", na/sum*100.0);
    println!("with R value: {}", best_r_value);
}