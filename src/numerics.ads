with Ada.Containers.Vectors, Interfaces.C, Interfaces.C.Pointers, Ada.Text_IO;
with Ada.Numerics.Generic_Elementary_Functions;
package Numerics is
   package C renames Interfaces.C;
   
   -------- Define types (Real, Int, Pos, Nat) ----------------------
   type Real is new C.double range C.double'First .. C.double'Last;
   type Int  is new C.long range C.long'First .. C.long'Last;
   subtype Pos is Int     range 0 .. Int'Last;
   subtype Nat is Pos     range 1 .. Pos'Last;
   
   -------- Define array types
   type Real_Array is array (Nat range <>) of aliased Real with Convention => C;
   type Int_Array  is array (Nat range <>) of aliased Int  with Convention => C;
   package Real_Functions is new Ada.Numerics.Generic_Elementary_Functions (Real);
   
   ------- Define Real_IO and Int_IO packages ------------------------
   package Int_IO is new Ada.Text_IO.Integer_IO (Int);
   package Real_IO is new Ada.Text_IO.Float_IO (Real);
   --  procedure New_Line (Spacing : in Ada.Text_IO.Positive_Count := 1) 
   --    renames Ada.Text_IO.New_Line;
   --  procedure Put_Line (Item : in String) renames Ada.Text_IO.Put_Line;
   --  procedure Put (Item : in String) renames Ada.Text_IO.Put;
   
   ------- Define Real_Vector and Int_Vector packages ------------------------
   package IV_Package is new Ada.Containers.Vectors (Nat, Int, "=");
   package RV_Package is new Ada.Containers.Vectors (Nat, Real, "=");
   subtype Int_Vector  is IV_Package.Vector;
   subtype Real_Vector is RV_Package.Vector;
   
   -- Vectorize & To_Array are needed in Triplet_To_Matrix
   function Vectorize (Item : in Real_Array) return Real_Vector;
   function Vectorize (Item : in Int_Array)  return Int_Vector;
   function To_Array (Item : in Real_Vector) return Real_Array;
   function To_Array (Item : in Int_Vector) return Int_Array;
   
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
   
   -------- Max and Abs_Max functions ------------------
   function Max_Int_Array (Item : in Int_Array) return Int with Inline => True;
   function Abs_Max_IA (Item : in Int_Array) return Int with Inline => True;
   function Max_Real_Array (Item : in Real_Array) return Real with Inline => True;
   function Abs_Max_RA (Item : in Real_Array) return Real with Inline => True;
   function Max (Item : in Int_Array) return Int renames Max_Int_Array;
   function Max (Item : in Real_Array) return Real renames Max_Real_Array;
   function Max (X : in Int_Vector) return Int;
   function Max (X : in Real_Vector) return Real;
   function Abs_Max (Item : in Int_Array) return Int renames Abs_Max_IA;
   function Abs_Max (Item : in Real_Array) return Real renames Abs_Max_RA;
   function Abs_Max (Item : in Real_Vector) return Real;
   
   
   ------- Dot Products ---------------------------------
   function Dot_Product (Left_I, Right_J : in Int_Array;
			 Left_X, Right_Y : in Real_Array) return Real;
   function Dot_Product_RV (X, Y : in Real_Vector) return Real;
   function Dot_Product (X, Y : in Real_Vector) return Real 
     renames Dot_Product_RV;
   
   
   
   -------- Binary Operators ---------------------------
   function Mult_Int_Array (Left, Right : in Int_Array) return Boolean;
   function Mult_R_RV (Left  : in Real;
		       Right : in Real_Vector) return Real_Vector;
   function Add_RV_RV (Left, Right : in Real_Vector) return Real_Vector;
   function Minus_RV_RV (Left, Right : in Real_Vector) return Real_Vector;
   function "*" (Left, Right : in Int_Array) return Boolean renames Mult_Int_Array;
   function "*" (Left  : in Real;
		 Right : in Real_Vector) return Real_Vector renames Mult_R_RV;
   function "-" (Left, Right : in Real_Vector) return Real_Vector renames Minus_RV_RV;
   function "+" (Left, Right : in Real_Vector) return Real_Vector renames Add_RV_RV;
   
   
end Numerics;
