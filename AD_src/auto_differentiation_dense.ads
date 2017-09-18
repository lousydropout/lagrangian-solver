with Numerics, Numerics.Dense_Matrices;
use  Numerics, Numerics.Dense_Matrices;
generic
   N : Nat;
package Auto_Differentiation_Dense is
   type Evaluation_Level is (Value, Gradient, Hessian);
   Level : Evaluation_Level := Hessian;
   
   type AD_Type is private;
   type AD_Vector is array (1 .. N) of AD_Type;
   
   function Var (X  : in Real;
		 I  : in Nat;
		 Dx : in Real := 1.0) return AD_Type with Pre => I <= N;
   
   function Const (X : in Real) return AD_Type;
   function Var (X : in Real_Vector) return AD_Vector
     with Pre => X'Length = N;
   
   function Val (X : in AD_Type) return Real;
   function Grad (X : in AD_Type) return Real_Vector;
   function Hessian (X : in AD_Type) return Real_Matrix;
   
   function "+" (X, Y : in AD_Type) return AD_Type;
   function "-" (X, Y : in AD_Type) return AD_Type;
   function "*" (X, Y : in AD_Type) return AD_Type;
   function "/" (X, Y : in AD_Type) return AD_Type;
   function "**" (X : in AD_Type; K : in Integer) return AD_Type;
   
   function Sin (X : in AD_Type) return AD_Type;
   function Cos (X : in AD_Type) return AD_Type;
   function Tan (X : in AD_Type) return AD_Type;
   function Exp (X : in AD_Type) return AD_Type;
   function Log (X : in AD_Type) return AD_Type; 
   
   function "+" (X : in AD_Type) return AD_Type is (X);
   function "-" (X : in AD_Type) return AD_Type;
   
   function "+" (X : in Real;
		 Y : in AD_Type) return AD_Type is (Const (X) + Y);
   function "-" (X : in Real;
		 Y : in AD_Type) return AD_Type is (Const (X) - Y);
   function "+" (X : in AD_Type;
		 Y : in Real) return AD_Type is (X + Const (Y));
   function "-" (X : in AD_Type;
		 Y : in Real) return AD_Type is (X - Const (Y));
   function "*" (X : in Real;
		 Y : in AD_Type) return AD_Type;
   function "*" (X : in AD_Type;
		 Y : in Real) return AD_Type is (Y * X);
   function "/" (X : in AD_Type;
		 Y : in Real) return AD_Type is ((1.0 / Y) * X)
     with Pre => Y /= 0.0;
   function "/" (X : in Real;
		 Y : in AD_Type) return AD_Type is (X * (Y ** (-1)));
   ------------- procedures ----------------------
   procedure Print (X : in AD_Type);
   
private
   
   type AD_Type is
      record
	 Val : Real;
	 Grad : Real_Vector (1 .. N);
	 Hessian : Real_Matrix (1 .. N, 1 .. N);
      end record;
   
   G0 : constant Real_Vector (1 .. N) := (others => 0.0);
   H0 : constant Real_Matrix (1 .. N, 1 .. N) := (others => (others => 0.0));
end Auto_Differentiation_Dense;
