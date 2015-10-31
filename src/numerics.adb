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
   begin
      X.Reserve_Capacity (To'Length);
      X.Set_Length (0);
      for Y of To loop
	 X.Append (Y);
      end loop;
   end Set;
   
   procedure Set (X  : in out Int_Vector;
		  To : in     Int_Array) is
   begin
      X.Reserve_Capacity (To'Length);
      X.Set_Length (0);
      for Y of To loop
	 X.Append (Y);
      end loop;
   end Set;
   
   function Vectorize (Item : in Real_Array) return Real_Vector is
      X : Real_Vector;
   begin
      X.Reserve_Capacity (Item'Length);
      X.Set_Length (0);
      for Y of Item loop
	 X.Append (Y);
      end loop;
      return X;
   end Vectorize;
   
   function Vectorize (Item : in Int_Array) return Int_Vector is
      X : Int_Vector;
   begin
      X.Reserve_Capacity (Item'Length);
      X.Set_Length (0);
      for Y of Item loop
	 X.Append (Y);
      end loop;
      return X;
   end Vectorize;
   
   function To_Array (Item : in Real_Vector) return Real_Array is
      Result : Real_Array (1 .. Nat (Item.Length));
      I : Nat := 1;
   begin
      for Y of Item loop
	 Result (I) := Y;
	 I := I + 1;
      end loop;
      return Result;
   end To_Array;
   
   function To_Array (Item : in Int_Vector) return Int_Array is
      Result : Int_Array (1 .. Nat (Item.Length));
      I : Nat := 1;
   begin
      for Y of Item loop
	 Result (I) := Y;
	 I := I + 1;
      end loop;
      return Result;
   end To_Array;
   
   
   
   ----- Vector and Array functions
   
   function Basis_Vector (I, N : in Int) return Real_Vector is
      use Ada.Containers;
      Result : Real_Vector := RV_Package.To_Vector (0.0, Count_Type (N));
   begin
      Result (I) := 1.0;
      return Result;
   end Basis_Vector;
   
   procedure Print (V : in Real_Vector) is
   begin
      for X of V loop
	 Real_IO.Put (X); Ada.Text_IO.New_Line;
      end loop;
   end Print;
   
   
   procedure Set_Length (V : in out Real_Vector;
			 N : in     Int) is
   begin
      V.Set_Length (Ada.Containers.Count_Type (N));
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
