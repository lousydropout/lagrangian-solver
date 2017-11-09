with Numerics, Chebyshev, Dense_AD, Dense_AD.Integrator;
use  Numerics, Chebyshev;
package Sb_Package is
   Convergence_Exception : exception;
   -----------------------------------------------
   N   : constant Nat  := 2;
   K   : constant Nat  := 12;
   α   : Real;
   -----------------------------------------------
   package E_Solver   is new Dense_AD (1); 
   package AD_Package is new Dense_AD (2 * N); 
   package Integrator is new AD_Package.Integrator (K);
   use AD_Package, Integrator;
   
   function Phi (R : in AD_Type) return AD_Type;
   function KE (Q : in AD_Vector) return AD_Type;
   function PE (Q : in AD_Vector) return AD_Type;
   
   function Lagrangian (T : in Real;
   			X : in Vector) return AD_Type;
   function Hamiltonian (T : in Real;
   			 X : in Vector) return AD_Type;
   function Get_IC (X : in Vector;
   		    E : in Real) return Vector;
   
   function X1 (T : in Real) return Real;
   function X1 (X : in Vector) return Real;
   function Y1 (T : in Real) return Real;
   function R13 (X : in Vector) return Real_Vector;
   function R13 (X : in Vector) return Real;
   function V13 (X : in Vector) return Real_Vector;
   function V13_New (X : in Vector) return Real_Vector;
   function New_Vel (X : in Vector) return Vector;
   function Sgn (X : in Real) return Real;
   
   function Find_State_At_Level (Level : in Real;
				 A : in Array_Of_Vectors;
				 T  : in Real;
				 Dt : in Real;
				 Lower : in out Real;
				 Upper : in out Real;
				 Func : not null access function (X : Vector)
				   return Real) return Variable;

end Sb_Package;
