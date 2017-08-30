with Numerics, Numerics.Sparse_Matrices;
use  Numerics, Numerics.Sparse_Matrices;

package Forward_AD is
   
   type AD_Type is tagged private;
   type AD_Vector is array (Nat range <>) of AD_Type;
   
   
   function Var (X    : in Real;
		 I, N : in Nat;
		 Dx   : in Real	:= 1.0) return AD_Type;
   function Zero (N : in Nat) return AD_Type;
   function Var (X	: in Real_Array;
		 Length	: in Nat;
		 Start	: in Nat := 1) return AD_Vector;
   function Var (X : in Real_Array) return AD_Vector is 
      (Var (X => X, Length => X'Length));
   
   
   function Val (X : in AD_Type) return Real;
   function Grad (X : in AD_Type) return Sparse_Vector;
   function Grad (X : in AD_Type) return Real_Array;
   function Length (X : in AD_Type) return Pos;
   function "+" (X, Y : in AD_Type) return AD_Type with Pre => X.Length = Y.Length;
   function "-" (X, Y : in AD_Type) return AD_Type with Pre => X.Length = Y.Length;
   function "*" (X, Y : in AD_Type) return AD_Type with Pre => X.Length = Y.Length;
   function "/" (X, Y : in AD_Type) return AD_Type with Pre => X.Length = Y.Length;
   function "**" (X : in AD_Type; N : in Pos) return AD_Type;
   function Sin (X : in AD_Type) return AD_Type;
   function Cos (X : in AD_Type) return AD_Type;
   function Tan (X : in AD_Type) return AD_Type;
   function Exp (X : in AD_Type) return AD_Type;
   function Log (X : in AD_Type) return AD_Type; 
   
   function "+" (X : in AD_Type) return AD_Type;
   function "-" (X : in AD_Type) return AD_Type;
   
   function "*" (Y : in Real; X : in AD_Type) return AD_Type;
   function "*" (X : in AD_Type; Y : in Real) return AD_Type is (Y * X);
   function "/" (X : in AD_Type; Y : in Real) return AD_Type with Pre => Y /= 0.0;
   
   function "*" (X : in Real; Y : in AD_Vector) return AD_Vector;
   function "*" (X : in Real_Matrix; Y : in AD_Vector) return AD_Vector;
   
   procedure Print (X : in AD_Type);
   
   
   -- functions for AD_Vectors
   function Jacobian_Matrix (X : in AD_Vector) return Sparse_Matrix;
   
   
private
   
   type AD_Type is tagged
      record
	 N    : Pos := 0;
	 Val  : Real;
	 Grad : Numerics.Sparse_Vector;
      end record;
   
end Forward_AD;
