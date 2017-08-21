with Numerics, Ada.Text_IO, Molecular_Dynamics;
use  Numerics, Ada.Text_IO, Molecular_Dynamics;

procedure Bubble is
   use Real_Functions, Real_IO, Int_IO;
   N  : constant Nat := 4;
   M  : constant Nat := 10;
   Dt : constant Real := 5.0e-5;
   U  : constant Pos2D := (1.0, 0.0);
   V  : constant Pos2D := (0.0, 1.0);
   T_Final : constant Real := 4.0;
   T_Output : constant Real := 0.1;
   Zero : constant Pos2D := (0.0, 0.0);
   
   File : File_Type;
   
   R, S, R_New : Pos2D_Vector (1 .. 2 * N * M + N + M) := (others => Zero);
   Is_BC : BC_Vectype (R'Range);
   
   Vac : constant Nat := 4 * (M + 1) + 4;
   T, Tmp : Real := 0.0;
   Iter : Nat := 1;
   
   
   
begin
   
   Put (Rand); New_Line;
   
   Initialize_Lattice (N, M, Vac, R, Is_BC, False);
   R_New := R;
   
   Create (File => File, Name => "out.xyz");
   
   -- Initializations ----------
   Output (R, File);

   for I in R'Range loop
      if Is_BC (I) = True then -- move particle I
   	 R_New (I) := R (I) + 0.1 * Dt * U * R(I).Y;
      end if;
   end loop;
   
   
   
   
   Output (R_New, File);
   -----------------------------
   while (T < T_Output) loop
      S := Verlet (R => R,
		   R_New => R_New,
		   Is_BC => Is_BC,
		   Dt => Dt);
      R     := R_New;
      R_New := S;
      T     := T + Dt;
   end loop;
   Output (S, File);
   
   
   while (T < T_Final) loop
      Tmp := Tmp + T_Output;
      
      while (T < Tmp) loop
	 S := Verlet (R => R,
		      R_New => R_New,
		      Is_BC => Is_BC,
		      Dt => Dt);
	 R     := R_New;
	 R_New := S;
	 T     := T + Dt;
      end loop;
      Put ("Iter = "); Put (Iter); New_Line; Iter := Iter + 1;
      Put ("Time: "); Put (T); New_Line;
      Output (S, File);
   end loop;
   
   
   
   Close (File => File);
   
end Bubble;
