with Ada.Numerics.Generic_Elementary_Functions, Ada.Numerics.Float_Random;
with Ada.Numerics, Ada.Containers.Vectors, Ada.Text_IO;
with Ada.Numerics.Generic_Real_Arrays;
package Numerics is
   π : constant := Ada.Numerics.π;
   
   
   -------- Define types (Real, Int, Pos, Nat) ----------------------
   type Real is new Long_Long_Float 
     range Long_Long_Float'First .. Long_Long_Float'Last;
   subtype Pos is Integer range 0 .. Integer'Last;
   subtype Nat is Integer range 1 .. Integer'Last;
   
   type Sparse_Vector is private;
   
   type Pos2D is record
      X, Y : Real;
   end record;
   
   type Pos2D_Vector is array (Nat range <>) of Pos2D;
   
   -------- Define array types -----------------------------------
   type Real_Vector is array (Nat range <>) of Real;
   type Int_Array  is array (Nat range <>) of Integer;
   type Real_Matrix is array (Nat range <>, Nat range <>) of Real;
   package Real_Functions is new Ada.Numerics.Generic_Elementary_Functions (Real);
   --  package RA is new Ada.Numerics.Generic_Real_Arrays (Real); 
   --  subtype Real_Matrix is RA.Real_Matrix;
   --  function Solve (A : in Real_Matrix;
   --  		   B : in Real_Vector) return Real_Vector
   --    with Pre => A'Length (1) = A'Length (2);
   
   ------- Define Real_IO and Int_IO packages ------------------------
   package Int_IO  is new Ada.Text_IO.Integer_IO (Integer);
   package Real_IO is new Ada.Text_IO.Float_IO   (Real);
   
   -------- Define random variable function -----------------------
   function Rand return Real;
   function Rand (N : in Nat) return Real_Vector;
   
   
   
   function Sparse (X	: in Real_Vector;
		    N	: in Pos	:= 0;
		    Tol	: in Real	:= 1.0e-20) return Sparse_Vector;
   
   -----------   Real_Vector functions ---------------------------
   function "-" (X : in Real_Vector) return Real_Vector;
   function "+" (X : in Real_Vector;
		 Y : in Real_Vector) return Real_Vector 
     with Pre => X'Length = Y'Length;
   function "-" (X : in Real_Vector;
		 Y : in Real_Vector) return Real_Vector
     with Pre => X'Length = Y'Length;
   
   function "*" (Left  : in Real;
		 Right : in Real_Vector) return Real_Vector;
   function "*" (Left  : in Real_Vector;
		 Right : in Real) return Real_Vector is (Right * Left);
   
   -----------   Real_Vector functions ---------------------------
   function "*" (A : in Real_Matrix;
		 X : in Real_Vector) return Real_Vector
     with Pre => A'Length (2) = X'Length;
   -------- Pos2D functions --------------------------------------
   function "+" (X : in Pos2D) return Pos2D is (X);
   function "-" (X : in Pos2D) return Pos2D;
   function "+" (X, Y : in Pos2D) return Pos2D;
   function "-" (X, Y : in Pos2D) return Pos2D;
   function Dot (X, Y : in Pos2D) return Real;
   function "*" (X, Y : in Pos2D) return Real renames Dot;

   function "*" (X : in Real;
		 Y : in Pos2D) return Pos2D is (X * Y.X, X * Y.Y);
   function "*" (X : in Pos2D;
		 Y : in Real) return Pos2D is (Y * X);
   
   function "+" (X, Y : in Pos2D_Vector) return Pos2D_Vector;
   function "-" (X, Y : in Pos2D_Vector) return Pos2D_Vector;
   function "-" (X : in Pos2D_Vector) return Pos2D_Vector;
   function "*" (X : in Real;
		 Y : in Pos2D_Vector) return Pos2D_Vector;
   function "*" (X : in Pos2D_Vector;
		 Y : in Real) return Pos2D_Vector is (Y * X);
   
   function Norm (X : in Pos2D) return Real;
   function To_Array (Xvec : in Pos2D_Vector) return Real_Vector;
   
   ------- Sparse_Vector Functions ----------------------------------------
   procedure Print (X : in Sparse_Vector); 
   function To_Array (X	  : in Sparse_Vector) return Real_Vector;
   procedure Set_Length (X : in out Sparse_Vector;
			 N : in     Pos);
   procedure Set (Item : in out Sparse_Vector;
   		  I    : in     Nat;
   		  X    : in     Real)
     with Pre => I <= Length (Item);
   procedure Add (Item : in out Sparse_Vector;
   		  I    : in     Nat;
   		  X    : in     Real)
     with Pre => I <= Length (Item);
   function "+" (A, B : in Sparse_Vector) return Sparse_Vector;
   function "*" (A : in Real;
		 B : in Sparse_Vector) return Sparse_Vector;
   function "*" (A : in Sparse_Vector;
		 B : in Real) return Sparse_Vector is (B * A);
   function "/" (A : in Sparse_Vector;
		 B : in Real) return Sparse_Vector is ((1.0 / B) * A);
   function "-" (A : in Sparse_Vector) return Sparse_Vector is ("*" (-1.0, A));
   function "-" (A : in Sparse_Vector) return Real_Vector is (To_Array (-A));
   function "-" (A, B : in Sparse_Vector) return Sparse_Vector is (A + (-B));
   
   
   function "*" (A : in Real_Matrix;
		 B : in Sparse_Vector) return Sparse_Vector;
   
   function Remove_1stN (X : in Sparse_Vector;
			 N : in Pos) return Sparse_Vector;
   
   ------- Norm --------------------------
   function Norm (X : in Sparse_Vector) return Real;
   function Length (X : in Sparse_Vector) return Pos;
   function Norm (X : in Real_Vector) return Real;
   
   -------- Max and Abs_Max functions ------------------
   function Max_Int_Array (Item : in Int_Array) return Integer;
   function Abs_Max_IA (Item : in Int_Array) return Integer;
   function Max_Real_Array (Item : in Real_Vector) return Real;
   function Abs_Max_RA (Item : in Real_Vector) return Real;
   
   function Max (Item : in Int_Array) return Integer renames Max_Int_Array;
   function Max (Item : in Real_Vector) return Real renames Max_Real_Array;
   function Abs_Max (Item : in Int_Array) return Integer renames Abs_Max_IA;
   function Abs_Max (Item : in Real_Vector) return Real renames Abs_Max_RA;
   
   
   ------- Dot Products ---------------------------------
   function Dot_Product (Left_I, Right_J : in Int_Array;
			 Left_X, Right_Y : in Real_Vector) return Real
     with Pre => Left_I'Length = Left_X'Length
     and Right_J'Length = Right_Y'Length;
   
   

   --  function "**" (Left : in Real; Right : in Int) return Real;
   function "or" (X, Y : in Sparse_Vector) return Sparse_Vector;
   
private
   
   Gen : Ada.Numerics.Float_Random.Generator;
   
   ------- Define RVector and IVector packages ------------------------
   package IV_Package is new Ada.Containers.Vectors (Nat, Integer,  "=");
   package RV_Package is new Ada.Containers.Vectors (Nat, Real, "=");
   
   subtype IVector  is IV_Package.Vector;
   subtype RVector is RV_Package.Vector;
   
   
   
   -- Vectorize & To_Array are needed in Triplet_To_Matrix
   procedure Set (X  : in out RVector;
		  To : in     Real_Vector);
   procedure Set (X  : in out IVector;
		  To : in     Int_Array);

   function Vectorize (Item : in Real_Vector) return RVector;
   function Vectorize (Item : in Int_Array)  return IVector;
   function To_Array (Item : in RVector) return Real_Vector;
   function To_Array (Item : in IVector) return Int_Array;

   ----- Vector and Array functions
   procedure Set_Length (V : in out RVector;
			 N : in     Integer);
   function Length (X : in RVector) return Integer;
   function Max (X : in IVector) return Integer;
   function Max (X : in RVector) return Real;
   type Sparse_Vector is record
      NMax : Pos := 0;
      X    : RVector;
      I    : IVector;
   end record;
   
end Numerics;
