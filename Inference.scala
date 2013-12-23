import breeze.linalg._

def anls1d(X: Vector[Int], W: Vector[Int], U: DenseMatrix[Double], V: DenseMatrix[Double], 
    alpha: Double = 0.05, Max_steps: Int=300, lambda1: Double = 0.015) : DenseVector[Double] = {

    val nzc = W.sum
    if (nzc == 0) return U, 0
    val MN = V.rows
    val GN = V.cols
    val X_reduced = for(i<-0 until MN if W[i] !=0) yield X[i].toDouble
    //val V_reduced = for(i<-0 until MN if W[i] !=0) yield X[i]
    V_reduced = DenseMatrix.zeros(nzc, GN)
    var j = 0
    for(i<-0 until MN if W[i]!=0){
        V_reduced(j, ::) := V(j, ::)
        j++
    }
    
    var u_new = U.copy
    var vv = V_reduced * V_reduced
    
    def RMSE1(X: Vector[Double], U: DenseMatrix[Double], V: DenseMatrix[Double]): Double = {
        var uv = U * V.t
        var temp = DenseMatrix(X.toArray) - U * V.t
        temp *= temp
        temp.data.sum.toDouble / temp.size
    }
    
    var lrmse = RMSE1(X_reduced, u_new, V_reduced)
    var u_best = u_new.copy
    var nzc1 = nzc.toDouble
    
    for(i<-0 to MAX_steps){
        lu = u_new.copy
        u_new += (X_reduced / nzc1 - vv/nzc1 - lambda1 * u_new) * alpha
        
    }
}