
object Hansen {
  
  /* A solvent with the three Hansen parameters. */
  case class Solvent(
    d_param: Double, 
    p_param: Double, 
    h_param: Double) 
  {
    def *(factor: Double) =
      Solvent(
        d_param * factor,
        p_param * factor,
        h_param * factor)
    
    def +(other: Solvent) =
      Solvent(
        d_param + other.d_param,
        p_param + other.p_param,
        h_param + other.h_param)
        
    def r_value(other: Solvent) =
      List(
        d_param - other.d_param,
        p_param - other.p_param,
        h_param - other.h_param)
        .map(x => x * x)
        .sum
  }
    
  object Solvent {
    
    def combine(proportionated: Seq[(Solvent, Double)]) = {
      val sumProp = proportionated.map(_._2).sum
      proportionated
        .map { case (s, p) => s * p }
        .reduce( (s1,s2) => s1 + s2 )
        .*(1/sumProp)
    }
    
    val ACETIC_ACID      = Solvent(6.8, 6.0, 9.2)
    
    val DIOCTYL_PHTALATE = Solvent(8.1, 3.4, 1.5)
    
    val NITROAMINE       = Solvent(8.6, 6.8, 0.0)
    
    /* Mystery solvent MEP */
    val MEP              = Solvent(7.8, 4.4, 2.5)
    
  }
  
  def main(args: Array[String]) {
    val prospective =
      Solvent.combine(List(
        (Solvent.ACETIC_ACID,      0.2),
        (Solvent.DIOCTYL_PHTALATE, 0.5),
        (Solvent.NITROAMINE,       0.3)
      ))
    
    print("R-value for 20% AA, 50% DP, 30% NA:")
    println(prospective.r_value(Solvent.MEP))
    
    var bestProportion = (0.0,0.0,0.0)
    var bestRvalue     = java.lang.Double.MAX_VALUE
    
    for (
        aa <- 0.0 to 1.0 by 0.01;
        dp <- 0.0 to 1.0 by 0.01;
        if (aa + dp <= 1.0);
        na  = 1.0 - aa - dp) {
      val prospective = 
        Solvent.combine(List(
          (Solvent.ACETIC_ACID,      aa),
          (Solvent.DIOCTYL_PHTALATE, dp),
          (Solvent.NITROAMINE,       na)))
      
      val currentR = prospective.r_value(Solvent.MEP)
      
      if (currentR < bestRvalue) {
        bestRvalue     = currentR
        bestProportion = (aa, dp, na)
      }
    }
    
    val (aa, dp, na) = bestProportion
    println("Best match for MEP is:")
    println(s"${aa*100}% Acetic Acid")
    println(s"${dp*100}% Dioctyl Phtalate")
    println(s"${na*100}% Nitroamine")
    println(s"with R value: ${bestRvalue}")
  }
  
}