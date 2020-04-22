SELECT breed_id, variant_id, breed_name 
FROM OMIA.Variant_Breed
INNER JOIN (
	SELECT * 
	FROM Breed 
	WHERE gb_species_id = 9796
    ) AS breeds 
USING(breed_id)