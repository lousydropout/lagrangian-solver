package body Molecular_Dynamics is
   
   function Verlet (R	  : in Pos2D_Vector;
		    R_New : in Pos2D_Vector;
		    Is_BC : in BC_Vectype;
		    Dt	  : in Real) return Pos2D_Vector is
      F : Pos2D_Vector := Calculate_Forces (R_New, Is_BC);
   begin
      return (2.0 * R_New - R + (Dt**2) * F);
   end Verlet;
   
   
   function LJ_Force (From : in Pos2D;
		      To   : in Pos2D) return Pos2D is
      use Real_Functions;
      E  : constant Real := 10000.0;
      R  : constant Real := Norm (To - From);
      R6 : constant Real := (1.0 / R) ** 6;
      F  : Real;
   begin
      F := (E / R ** 2) * R6 * (R6 - 1.0);
      return (F * (To - From));
   end LJ_Force;
   
   
   
   function Calculate_Forces (R	    : in Pos2D_Vector;
			      Is_BC : in BC_Vectype) return Pos2D_Vector is
      Force : Pos2D := (0.0, 0.0);
      F     : Pos2D_Vector (R'Range) := (others => (0.0, 0.0));
   begin
      for I in R'First .. R'Last loop
	 for J in I + 1 .. R'Last loop
	    Force := LJ_Force (From => R (I), To => R (J));
	    F (I) := F (I) - Force;
	    F (J) := F (J) + Force;
	 end loop;
      end loop;
      
      Force := (0.0, 0.0);
      for X of F loop
	 Force := Force + X;
      end loop;
      if Norm (Force) > 0.1 then
	 Put (Norm (Force));
	 Put_Line ("Error");
      end if;
      
      for I in R'Range loop
	 if Is_BC (I) = True then F (I) := (0.0, 0.0); end if;
      end loop;
      
      return F;
   end Calculate_Forces;
   
   
   
   
   
   
   procedure Initialize_Lattice (N, M	 : in     Nat;
				 Vac	 : in     Nat;
				 Lattice :    out Pos2d_Vector;
				 Is_BC	 :    out BC_Vectype;
				 Sides	 : in     Boolean      := False) is
      K, Ind : Nat := 1;
      Ratio : Real := 1.0;
   begin
      Is_BC := (others => False);
      Is_BC (1 .. 2 * M + 2) := (others => True);
      
      -- Row 2N
      for J in 0 .. M loop
	 Lattice (K).X := Ratio * Real (J);
	 Lattice (K).Y := Ratio * Real (2 * N) * Sqrt (3.0) * 0.5;
	 K := K + 1; Ind := Ind + 1;
      end loop;
      
      -- Row 0 
      for J in 0 .. M loop
	 Lattice (K).X := Ratio * Real (J);
	 Lattice (K).Y := 0.0;
	 K := K + 1; Ind := Ind + 1;
      end loop;
      
      -- The in-between rows
      for I in 1 .. 2 * N - 1 loop
	 if Sides = True then Is_BC (K) := True; end if;
	    
	 if I mod 2 = 0 then -- if row I is even
	    for J in 0 .. M loop
	       if Ind /= Vac then
		  Lattice (K).X := Ratio * Real (J);
		  Lattice (K).Y := Ratio * Real (I) * Sqrt (3.0) * 0.5;
		  K := K + 1; 
	       end if;
	       Ind := Ind + 1;
	    end loop;
	 else -- if row I is odd
	    for J in 1 .. M loop
	       if Ind /= Vac then
	    	  Lattice (K).X := Ratio * Real (J) - 0.5;
	    	  Lattice (K).Y := Ratio * Real (I) * Sqrt (3.0) * 0.5;
	    	  K := K + 1; 
	       end if;
	       Ind := Ind + 1;
	    end loop;
	    null;
	 end if;
	 if Sides = True then Is_BC (K - 1) := True; end if;
	 
      end loop;
      
	 for I in Lattice'Range loop
	    --  if Is_BC (I) = False then
	       Lattice (I).X := Lattice (I).X; --  + 0.05 * Rand;
	       Lattice (I).Y := Lattice (I).Y; -- + 0.05 * Rand;
	    --  end if;
	 end loop;

   end Initialize_Lattice;
      

   
   
   
   
   
   
   
   procedure Output (R	  : in Pos2D_Vector;
		     File : in File_Type) is
   begin
      Put (File, R'Length); New_Line (File);
      Put_Line (File, "Properties=pos:R:2");
      
      for Item of R loop
	 Put (File => File, Item => Item.X);
	 Put (File => File, Item => "     ");
	 Put (File => File, Item => Item.Y);
	 Put (File => File, Item => "     ");
	 Put (File => File, Item => "0.5");
	 New_Line (File => File);
      end loop;
   end Output;
   

   
end Molecular_Dynamics;
