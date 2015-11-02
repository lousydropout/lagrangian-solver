with Ada.Containers.Vectors, Interfaces.C, Interfaces.C.Pointers, Ada.Text_IO;
with Ada.Numerics.Generic_Elementary_Functions, Ada.Numerics.Float_Random;
with Ada.Text_IO;
package Numerics is
   package C renames Interfaces.C;
   
   -------- Define types (Real, Int, Pos, Nat) ----------------------
   type Real is new C.double range C.double'First .. C.double'Last;
   type Int  is new C.long range C.long'First .. C.long'Last;
   subtype Pos is Int     range 0 .. Int'Last;
   subtype Nat is Int     range 1 .. Int'Last;
   
   type Sparse_Vector is tagged private;
   
   -------- Define random variable function -----------------------
   function Rand return Real;
   
   -------- Define array types -----------------------------------
   type Real_Array is array (Nat range <>) of aliased Real with Convention => C;
   type Int_Array  is array (Nat range <>) of aliased Int  with Convention => C;
   type Real_Matrix is array (Nat range <>, Nat range <>) of aliased Real with Convention => C;
   package Real_Functions is new Ada.Numerics.Generic_Elementary_Functions (Real);
   
   ------- Define Real_IO and Int_IO packages ------------------------
   package Int_IO is new Ada.Text_IO.Integer_IO (Int);
   package Real_IO is new Ada.Text_IO.Float_IO (Real);
   --  procedure New_Line (Spacing : in Ada.Text_IO.Positive_Count := 1) 
   --    renames Ada.Text_IO.New_Line;
   --  procedure Put_Line (Item : in String) renames Ada.Text_IO.Put_Line;
   --  procedure Put (Item : in String) renames Ada.Text_IO.Put;
   
   ------- Define Real_Vector and Int_Vector packages ------------------------
   package IV_Package is new Ada.Containers.Vectors (Nat, Int,  "=");
   package RV_Package is new Ada.Containers.Vectors (Nat, Real, "=");
   subtype Int_Vector  is IV_Package.Vector;
   subtype Real_Vector is RV_Package.Vector;
   
   -- Vectorize & To_Array are needed in Triplet_To_Matrix
   procedure Set (X  : in out Real_Vector;
		  To : in     Real_Array);
   procedure Set (X  : in out Int_Vector;
		  To : in     Int_Array);

   function Vectorize (Item : in Real_Array) return Real_Vector;
   function Vectorize (Item : in Int_Array)  return Int_Vector;
   function To_Array (Item : in Real_Vector) return Real_Array;
   function To_Array (Item : in Int_Vector) return Int_Array;
   function Sparse (X   : in Real_Vector;
		    N   : in Pos	:= 0;
		    Tol	: in Real	:= 1.0e-20) return Sparse_Vector;
   function Sparse (X	: in Real_Array;
		    N	: in Pos	:= 0;
		    Tol	: in Real	:= 1.0e-20) return Sparse_Vector;
   
   ------- Sparse_Vector Functions ----------------------------------------
   procedure Print (X : in Sparse_Vector); 
   function To_Array (X	  : in Sparse_Vector) return Real_Array;
   procedure Set_Length (X : in out Sparse_Vector;
			 N : in     Pos);
   procedure Set (Item : in out Sparse_Vector;
   		  I    : in     Nat;
   		  X    : in     Real);
   function "+" (A, B : in Sparse_Vector) return Sparse_Vector;
   function "*" (A : in Real;
		 B : in Sparse_Vector) return Sparse_Vector;
   function "*" (A : in Sparse_Vector;
		 B : in Real) return Sparse_Vector is (B * A);
   function "/" (A : in Sparse_Vector;
		 B : in Real) return Sparse_Vector is ((1.0 / B) * A);
   function "-" (A : in Sparse_Vector) return Sparse_Vector is ("*"(-1.0, A));
   function "-" (A, B : in Sparse_Vector) return Sparse_Vector is (A + (-B));
   
   
   ----- Vector and Array functions
   function Basis_Vector (I, N : in Int) return Real_Vector;
   procedure Print (V : in Real_Vector);
   procedure Set_Length (V : in out Real_Vector;
			 N : in     Int);
   function Length (X : in Real_Vector) return Int;
   
   ------- Norm --------------------------
   function Norm2_RV (X : in Real_Vector) return Real;
   function Norm_RV (X : in Real_Vector) return Real;
   function Norm2 (X : in Real_Vector) return Real renames Norm2_RV;
   function Norm (X : in Real_Vector) return Real renames Norm_RV;
   function Norm (X : in Sparse_Vector) return Real;
   function Length (X : in Sparse_Vector) return Pos;
   function Norm (X : in Real_Array) return Real;
   
   -------- Max and Abs_Max functions ------------------
   function Max_Int_Array (Item : in Int_Array) return Int;
   function Abs_Max_IA (Item : in Int_Array) return Int;
   function Max_Real_Array (Item : in Real_Array) return Real;
   function Abs_Max_RA (Item : in Real_Array) return Real;
   function Max (Item : in Int_Array) return Int renames Max_Int_Array;
   function Max (Item : in Real_Array) return Real renames Max_Real_Array;
   function Max (X : in Int_Vector) return Int;
   function Max (X : in Real_Vector) return Real;
   function Abs_Max (Item : in Int_Array) return Int renames Abs_Max_IA;
   function Abs_Max (Item : in Real_Array) return Real renames Abs_Max_RA;
   function Abs_Max (Item : in Real_Vector) return Real;
   
   
   ------- Dot Products ---------------------------------
   function Dot_Product (Left_I, Right_J : in Int_Array;
			 Left_X, Right_Y : in Real_Array) return Real
     with Pre => Left_I'Length = Left_X'Length
     and Right_J'Length = Right_Y'Length;
   function Dot_Product_RV (X, Y : in Real_Vector) return Real
     with Pre => Pos (X.Length) = Pos (Y.Length);
   function Dot_Product (X, Y : in Real_Vector) return Real 
     renames Dot_Product_RV;
   
   
   
   -------- Binary Operators ---------------------------
   --  function Mult_Int_Array (Left, Right : in Int_Array) return Boolean;
   function Mult_R_RV (Left  : in Real;
		       Right : in Real_Vector) return Real_Vector;
   function Add_RV_RV (Left, Right : in Real_Vector) return Real_Vector
     with Pre => Pos (Left.Length) = Pos (Right.Length);
   function Minus_RV_RV (Left, Right : in Real_Vector) return Real_Vector
     with Pre => Pos (Left.Length) = Pos (Right.Length);
   function "*" (Left  : in Real;
		 Right : in Real_Vector) return Real_Vector renames Mult_R_RV;
   function "-" (Left, Right : in Real_Vector) return Real_Vector renames Minus_RV_RV;
   function "+" (Left, Right : in Real_Vector) return Real_Vector renames Add_RV_RV;
   function "*" (A : in Real_Vector;
		 B : in Real) return Real_Vector is (B * A);
   function "/" (A : in Real_Vector;
		 B : in Real) return Real_Vector is ((1.0 / B) * A);
   
private
   
   Gen : Ada.Numerics.Float_Random.Generator;
   
   type Sparse_Vector is tagged record
      NMax : Pos := 0;
      X    : Real_Vector;
      I    : Int_Vector;
   end record;
   
end Numerics;
