SELECT *, 
	CASE 
		WHEN LOWER(reference_sequence) LIKE "equcab%" THEN "horse"
        WHEN conf_horse = "horse" THEN "horse"
        ELSE null
	END AS horse_filter
		
FROM OMIA.Variant
LEFT JOIN (

	SELECT variant_id, GROUP_CONCAT(breed_id) AS breed_id, GROUP_CONCAT(breed_name) AS breed_name, "horse" AS conf_horse
	FROM OMIA.Variant_Breed
	INNER JOIN (
		SELECT * 
		FROM Breed 
		WHERE gb_species_id = 9796
		) AS breeds 
	USING(breed_id)
    GROUP BY variant_id
    
) AS breed_variants USING(variant_id)

LEFT JOIN (
	SELECT gene_id, GROUP_CONCAT(synonym) AS gene
	FROM OMIA.GeneSynonym
    GROUP BY gene_id
    # gene_ids may appears twice due to gene name synonyms
    ) AS genes USING(gene_id)

LEFT JOIN(
	SELECT variant_type_id as type_id, variant_type_name 
	FROM OMIA.Variant_Type
) AS var_types USING(type_id)

HAVING horse_filter = "horse"


