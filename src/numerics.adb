with Ada.Numerics.Generic_Real_Arrays;
package body Numerics is
   
   package Real_Arrays is new Ada.Numerics.Generic_Real_Arrays (Real);
   
   function Rand return Real is 
      use Ada.Numerics.Float_Random;
   begin
      return Real (Random (Gen));
   end Rand;
   
   -- Vectorize & To_Array are needed in Triplet_To_Matrix
   procedure Set (X  : in out Real_Vector;
		  To : in     Real_Array) is
      Offset : constant Int := To'First - 1;
   begin
      X.Set_Length (To'Length);
      for K in 1 .. Int (To'Length) loop
	 X (K) := To (K + Offset);
      end loop;
   end Set;
   
   procedure Set (X  : in out Int_Vector;
		  To : in     Int_Array) is
      Offset : constant Int := To'First - 1;
   begin
      X.Set_Length (To'Length);
      for K in 1 .. Int (To'Length) loop
	 X (K) := To (K + Offset);
      end loop;
   end Set;
   
   function Vectorize (Item : in Real_Array) return Real_Vector is
      Vector : Real_Vector;
      Offset : constant Int := Item'First - 1;
   begin
      Vector.Set_Length (Item'Length);
      for K in 1 .. Int (Item'Length) loop
   	 Vector (K) := Item (K + Offset);
      end loop;
      return Vector;
   end Vectorize;
   
   function Vectorize (Item : in Int_Array) return Int_Vector is
      Vector : Int_Vector;
      Offset : constant Int := Item'First - 1;
   begin
      Vector.Set_Length (Item'Length);
      for K in 1 .. Int (Item'Length) loop
   	 Vector (K) := Item (K + Offset);
      end loop;
      return Vector;
   end Vectorize;
   
   function To_Array (Item : in Real_Vector) return Real_Array is
      Result : Real_Array (1 .. Nat (Item.Length));
   begin
      for K in Result'Range loop
	 Result (K) := Item (K);
      end loop;
      return Result;
   end To_Array;
   
   function To_Array (Item : in Int_Vector) return Int_Array is
      Result : Int_Array (1 .. Nat (Item.Length));
   begin
      for K in Result'Range loop
	 Result (K) := Item (K);
      end loop;
      return Result;
   end To_Array;
   
   
   
   ----- Vector and Array functions
   
   function Basis_Vector (I, N : in Int) return Real_Vector is
      use Ada.Containers;
      Result : Real_Vector;
   begin
      Result.Set_Length (Count_Type (N));
      Result (I) := 1.0;
      return Result;
   end Basis_Vector;
   
   procedure Print (V : in Real_Vector) is
      use Real_IO, Ada.Text_IO;
   begin
      for X of V loop
	 Put (X); New_Line;
      end loop;
   end Print;
   
   
   procedure Set_Length (V : in out Real_Vector;
			 N : in     Int) is
      use Ada.Containers;
   begin
      V.Set_Length (Count_Type (N));
   end Set_Length;

   function Length (X : in Real_Vector) return Int is (Int (X.Length));
   
   
   
   
   ------- Norm --------------------------
   
   function Norm2_RV (X : in Real_Vector) return Real is separate;
   function Norm_RV (X : in Real_Vector) return Real is separate;
   
   
   
   -------- Max and Abs_Max functions ------------------
   
   function Max_Int_Array (Item : in Int_Array) return Int is separate;
   function Max_Real_Array (Item : in Real_Array) return Real is separate;
   function Abs_Max_IA (Item : in Int_Array) return Int is separate;
   function Abs_Max_RA (Item : in Real_Array) return Real is separate;
   
   function Max (X : in Int_Vector) return Int is
      Result : Int;
   begin
      case X.Length is
	 when 0 => return 0;
	 when 1 => return X (1);
	 when others =>
	    Result := X (1);
	    for I in 2 .. Int (X.Length) loop
	       if X (I) > Result then Result := X (I); end if;
	    end loop;
	    return Result;
      end case;
   end Max;
   
   function Max (X : in Real_Vector) return Real is
      Result : Real;
   begin
      case X.Length is
	 when 0 => return 0.0;
	 when 1 => return X (1);
	 when others =>
	    Result := X (1);
	    for I in 2 .. Int (X.Length) loop
	       if X (I) > Result then Result := X (I); end if;
	    end loop;
	    return Result;
      end case;
   end Max;


   function Abs_Max (Item : in Real_Vector) return Real is
      Result : Real := 0.0;
   begin
      for X of Item loop
	 if abs (X) > Result then Result := X; end if;
      end loop;
      return Result;
   end Abs_Max;
   
   
   
   -------- Binary Operators ---------------------------
   function Dot_Product (Left_I, Right_J : in Int_Array;
			 Left_X, Right_Y : in Real_Array) return Real is separate;
   function Dot_Product_RV (X, Y : in Real_Vector) return Real is separate;

   function Mult_Int_Array (Left, Right : in Int_Array) return Boolean is separate;
   function Mult_R_RV (Left  : in Real;
		       Right : in Real_Vector) return Real_Vector is separate;
   function Add_RV_RV (Left, Right : in Real_Vector) return Real_Vector is separate;
   function Minus_RV_RV (Left, Right : in Real_Vector) return Real_Vector is separate;
   

begin
   Ada.Numerics.Float_Random.Reset (Gen);
end Numerics;
