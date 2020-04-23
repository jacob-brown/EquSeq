SELECT *, 
	CASE 
		WHEN LOWER(reference_sequence) LIKE "equcab%" THEN "horse"
        WHEN conf_horse = "horse" THEN "horse"
        ELSE null
	END AS horse_filter
		
FROM OMIA.Variant
LEFT JOIN (

	SELECT breed_id, variant_id, breed_name, "horse" AS conf_horse
	FROM OMIA.Variant_Breed
	INNER JOIN (
		SELECT * 
		FROM Breed 
		WHERE gb_species_id = 9796
		) AS breeds 
	USING(breed_id)
) AS breed_variants USING(variant_id)

LEFT JOIN (
	SELECT gene_id, synonym AS gene
	FROM OMIA.GeneSynonym
    ) AS genes USING(gene_id)

LEFT JOIN(
	SELECT variant_type_id as type_id, variant_type_name
	FROM OMIA.Variant_Type
) AS var_types USING(type_id)

HAVING horse_filter = "horse"




