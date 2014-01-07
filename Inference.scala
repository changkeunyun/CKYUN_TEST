import breeze.linalg._
import math._
//import scala.util.control.Breaks._

def RMSE1(X: DenseMatrix[Double], U: DenseMatrix[Double], V: DenseMatrix[Double]): Double = {
    var temp = X - U * V.t
    temp *= temp
    sqrt(temp.data.sum.toDouble / temp.size)
}

def anls1d(X: DenseVector[Int], W: DenseVector[Int], U: DenseMatrix[Double], V: DenseMatrix[Double], 
    alpha: Double = 0.05, MAX_steps: Int=300, lambda1: Double = 0.015) : DenseMatrix[Double] = {

    val t = System.nanoTime

    val nzc = W.sum.toInt
    //if (nzc == 0) return (U, 0.0)
    if (nzc == 0) return U
    val MN = V.rows
    val GN = V.cols
    var X_reduced_Array = for(i<-0 until MN if W(i) !=0) yield X(i).toDouble
    var X_reduced = DenseMatrix(X_reduced_Array.toArray)
    //val V_reduced = for(i<-0 until MN if W[i] !=0) yield X[i]
    var V_reduced = DenseMatrix.zeros[Double](nzc, GN)
    var j = 0
    for(i<-0 until MN if W(i)!=0){
        V_reduced(j, ::) := V(j, ::)
        j += 1
    }

    var u_new = U.copy
    var vv = V_reduced.t * V_reduced

    var lrmse = RMSE1(X_reduced, u_new, V_reduced)
    var brmse = lrmse
    var u_best = u_new.copy
    var nzc1 = nzc.toDouble

    var alpha1 = alpha
    
    breakable{
        for(i<-0 to MAX_steps){
            var lu = u_new.copy
            u_new += (X_reduced*V_reduced / nzc1 - u_new*vv/nzc1 - u_new * lambda1) * alpha1

            for(j <-0 until GN if u_new(0,j) < 0.0) u_new(0,j) = 0.0
            var nrmse = RMSE1(X_reduced, u_new, V_reduced)

            if(nrmse < brmse){
                brmse = nrmse
                u_best = u_new.copy
                //println(brmse)
            }
        
            if(lrmse > nrmse){
                alpha1 = alpha1 * 1.05
                if(lrmse - nrmse < 0.00001)
                    break()
                lrmse = nrmse
            }
            else{
                u_new = lu.copy
                alpha1 = alpha1 * 0.5
                if(alpha1 < 1e-10) 
                    break()            
            }
        }
    }
    println(System.nanoTime - t)
    //(u_new, brmse)
    println(brmse)
    u_new
}