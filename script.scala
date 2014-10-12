
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
    
    val aceticAcid = Solvent(6.8, 6.0, 9.2)
    
    val dioctylPhtalate = Solvent(8.1, 3.4, 1.5)
    
    val nitroamine = Solvent(8.6, 6.8, 0.0)
    
    /* Mystery solvent MEP */
    val mep = Solvent(7.8, 4.4, 2.5)
    
  }
  
  def main(args: Array[String]) {
    val prospective =
      Solvent.combine(List(
        (Solvent.aceticAcid,      0.2),
        (Solvent.dioctylPhtalate, 0.5),
        (Solvent.nitroamine,      0.3)
      ))
    
    println(s"R-value for 20% AA, 50% DP, 30% NA: ${
        prospective.r_value(Solvent.mep)}")
    
    val props = 0.0 to 1.0 by 0.01
    var bestProportion = (0.0,0.0,0.0)
    var bestRvalue = java.lang.Double.MAX_VALUE
    
    for (
        aa <- props;
        dp <- props;
        if (aa + dp <= 1.0);
        na  = 1.0 - aa - dp) {
      val prospective = 
        Solvent.combine(List(
          (Solvent.aceticAcid,      aa),
          (Solvent.dioctylPhtalate, dp),
          (Solvent.nitroamine,      na)))
      
      val currentR = prospective.r_value(Solvent.mep)
      
      if (currentR < bestRvalue) {
        bestRvalue = currentR
        bestProportion = (aa, dp, na)
      }
    }
    
    val (aa, dp, na) = bestProportion
    println(s"Best match for MEP is ${aa*100}% AA, ${dp*100
        }% DP and ${na*100}% NA with r_value: ${bestRvalue}")
  }
  
}