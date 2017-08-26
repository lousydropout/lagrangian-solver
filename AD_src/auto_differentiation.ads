with Numerics, Numerics.Sparse_Matrices;
use  Numerics, Numerics.Sparse_Matrices;

package Auto_Differentiation is
   type Evaluation_Level is (Value, Gradient, Hessian);
   Level : Evaluation_Level := Hessian;
   
   
   type AD_Type is private;
   type AD_2D is array (1 .. 2) of AD_Type;
      
   type AD_Vector is array (Nat range <>) of AD_Type;
   
   function Var (X    : in Real;
   		 I, N : in Nat;
   		 Dx   : in Real	:= 1.0) return AD_Type;
   function Const (X : in Real;
		   N : in Nat) return AD_Type;
   function Zero (N : in Nat) return AD_Type;
   function Var (X	: in Real_Array;
   		 Length	: in Nat;
   		 Start	: in Nat := 1) return AD_Vector;
   function Var (X : in Real_Array) return AD_Vector is 
      (Var (X => X, Length => X'Length));
   
   function Val (X : in AD_Type) return Real;
   function Grad (X : in AD_Type) return Sparse_Vector;
   function Grad (X : in AD_Type) return Real_Array;
   function Hessian (X : in AD_Type) return Sparse_Matrix;
   function Length (X : in AD_Type) return Pos;


   function "+" (X, Y : in AD_Type) return AD_Type;
   function "-" (X, Y : in AD_Type) return AD_Type;
   function "*" (X, Y : in AD_Type) return AD_Type;
   function "/" (X, Y : in AD_Type) return AD_Type;
   function "**" (X : in AD_Type; N : in Pos) return AD_Type;
   
   function Sin (X : in AD_Type) return AD_Type;
   function Cos (X : in AD_Type) return AD_Type;
   function Tan (X : in AD_Type) return AD_Type;
   function Exp (X : in AD_Type) return AD_Type;
   function Log (X : in AD_Type) return AD_Type; 
   
   function "+" (X : in AD_Type) return AD_Type is (X);
   function "-" (X : in AD_Type) return AD_Type;
   
   function "*" (Y : in Real; X : in AD_Type) return AD_Type;
   function "*" (X : in AD_Type; Y : in Real) return AD_Type is (Y * X);
   function "/" (X : in AD_Type; Y : in Real) return AD_Type  is ((1.0 / Y) * X)
     with Pre => Y /= 0.0;
   
   function "+" (X, Y : in AD_2D) return AD_2D;
   function "-" (X, Y : in AD_2D) return AD_2D;
   function Dot (X, Y : in AD_2D) return AD_Type;
   ------------- procedures ----------------------
   procedure Print (X : in AD_Type);
     
private
   
   type AD_Type is
      record
	 N       : Pos := 0;
	 Val     : Real;
	 Grad    : Sparse_Vector;
	 Hessian : Sparse_Matrix;
      end record;
   
   G0 : constant Sparse_Vector := Sparse ((0.0, 0.0));
   H0 : constant Sparse_Matrix := Zero (1);
end Auto_Differentiation;

