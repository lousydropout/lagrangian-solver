with Ada.Containers.Vectors, Interfaces.C, Interfaces.C.Pointers, Ada.Text_IO;
with Ada.Numerics.Generic_Elementary_Functions, Ada.Numerics.Float_Random;
with Ada.Numerics, Ada.Text_IO;
package Numerics is
   π : constant := Ada.Numerics.π;
   
   package C renames Interfaces.C;
   
   
   -------- Define types (Real, Int, Pos, Nat) ----------------------
   type Real is new C.double range C.double'First .. C.double'Last;
   type Int  is new C.long   range C.long'First   .. C.long'Last;
   subtype Pos is Int        range 0              .. Int'Last;
   subtype Nat is Int        range 1              .. Int'Last;
   
   ------- Define Real_Vector and Int_Vector packages ------------------------
   package IV_Package is new Ada.Containers.Vectors (Nat, Int,  "=");
   package RV_Package is new Ada.Containers.Vectors (Nat, Real, "=");
   subtype Int_Vector  is IV_Package.Vector;
   subtype Real_Vector is RV_Package.Vector;
   
   type Sparse_Vector is private;
   
   type Pos2D is record
      X, Y : Real;
   end record;
   
   type Pos2D_Vector is array (Nat range <>) of Pos2D;
   
   -------- Define array types -----------------------------------
   type Real_Array is array (Nat range <>) of aliased Real with Convention => C;
   type Int_Array  is array (Nat range <>) of aliased Int  with Convention => C;
   type Real_Matrix is array (Nat range <>, Nat range <>) of aliased Real with Convention => C;
   package Real_Functions is new Ada.Numerics.Generic_Elementary_Functions (Real);
   
   ------- Define Real_IO and Int_IO packages ------------------------
   package Int_IO  is new Ada.Text_IO.Integer_IO (Int);
   package Real_IO is new Ada.Text_IO.Float_IO   (Real);
   
   -------- Define random variable function -----------------------
   function Rand return Real;
   function Rand (N : in Nat) return Real_Array;
   
   
   
   function Sparse (X	: in Real_Array;
		    N	: in Pos	:= 0;
		    Tol	: in Real	:= 1.0e-20) return Sparse_Vector;
   
   -----------   Real_Array functions ---------------------------
   function "-" (X : in Real_Array) return Real_Array;
   function "+" (X : in Real_Array;
		 Y : in Real_Array) return Real_Array 
     with Pre => X'Length = Y'Length;
   function "-" (X : in Real_Array;
		 Y : in Real_Array) return Real_Array
     with Pre => X'Length = Y'Length;
   
   function "*" (Left  : in Real;
		 Right : in Real_Array) return Real_Array;
   function "*" (Left  : in Real_Array;
		 Right : in Real) return Real_Array is (Right * Left);
   
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
   function To_Array (Xvec : in Pos2D_Vector) return Real_Array;
   
   ------- Sparse_Vector Functions ----------------------------------------
   procedure Print (X : in Sparse_Vector); 
   function To_Array (X	  : in Sparse_Vector) return Real_Array;
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
   function "-" (A : in Sparse_Vector) return Real_Array is (To_Array (-A));
   function "-" (A, B : in Sparse_Vector) return Sparse_Vector is (A + (-B));
   
   
   function "*" (A : in Real_Matrix;
		 B : in Sparse_Vector) return Sparse_Vector;
   

   
   ------- Norm --------------------------
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
   function Abs_Max (Item : in Int_Array) return Int renames Abs_Max_IA;
   function Abs_Max (Item : in Real_Array) return Real renames Abs_Max_RA;
   
   
   ------- Dot Products ---------------------------------
   function Dot_Product (Left_I, Right_J : in Int_Array;
			 Left_X, Right_Y : in Real_Array) return Real
     with Pre => Left_I'Length = Left_X'Length
     and Right_J'Length = Right_Y'Length;
   
   

   --  function "**" (Left : in Real; Right : in Int) return Real;

   
private
   
   Gen : Ada.Numerics.Float_Random.Generator;
   
   
   
   
   -- Vectorize & To_Array are needed in Triplet_To_Matrix
   procedure Set (X  : in out Real_Vector;
		  To : in     Real_Array);
   procedure Set (X  : in out Int_Vector;
		  To : in     Int_Array);

   function Vectorize (Item : in Real_Array) return Real_Vector;
   function Vectorize (Item : in Int_Array)  return Int_Vector;
   function To_Array (Item : in Real_Vector) return Real_Array;
   function To_Array (Item : in Int_Vector) return Int_Array;

   ----- Vector and Array functions
   procedure Set_Length (V : in out Real_Vector;
			 N : in     Int);
   function Length (X : in Real_Vector) return Int;
   function Max (X : in Int_Vector) return Int;
   function Max (X : in Real_Vector) return Real;
   
   type Sparse_Vector is record
      NMax : Pos := 0;
      X    : Real_Vector;
      I    : Int_Vector;
   end record;
   
end Numerics;
