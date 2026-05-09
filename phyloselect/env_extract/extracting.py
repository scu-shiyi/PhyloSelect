import pandas as pd
import rasterio
from pathlib import Path
import geopandas as gpd
import os

def extract_raster_values(file_path, coords):
	with rasterio.open(file_path) as src:
		return [val[0] for val in src.sample(coords)]

def extracting_envs(lon_lat_data, env_type,database_dir):
	df = pd.read_csv(lon_lat_data, sep=',')
	df.columns = (df.columns.astype(str)
				  .str.replace("\ufeff", "", regex=False)
				  .str.strip()
				  .str.lower())
	if "lon" not in df.columns or "lat" not in df.columns:
		if "longitude" in df.columns and "latitude" in df.columns:
			df.rename(columns={"longitude": "lon", "latitude": "lat"}, inplace=True)
		else:
			raise ValueError(f"表格里必须包含 lat 和 lon 列！当前列名：{list(df.columns)}")

	coords = [(x, y) for x, y in zip(df["lon"], df["lat"])]
	env_type_map = {
		"bio19": "wc2.1_30s_bio",
		"ele": "wc2.1_30s_elev",
		"tavg": "wc2.1_30s_tavg",
		"wind": "wc2.1_30s_wind",
		"vapr": "wc2.1_30s_vapr",
		"prec": "wc2.1_30s_prec",
		"srad": "wc2.1_30s_srad"
	}
	if env_type.lower() not in env_type_map:
		raise ValueError(f"Unsupported environment type: {env_type}")
	data_dir = os.path.join(database_dir, env_type_map[env_type.lower()])

	if env_type.lower() == "bio19":
		for i in range(1, 20):
			bio_tif = os.path.join(data_dir, f"wc2.1_30s_bio_{i}.tif")
			if not os.path.exists(bio_tif):
				print(f"[Warning] Could not find bio{i} data")
				continue
			values = extract_raster_values(bio_tif, coords)
			df[f"bio_{i}"] = values
			print(f"[Info] bio{i} data extracted")

	elif env_type.lower() == "ele":
		ele_tif = os.path.join(data_dir, f"wc2.1_30s_elev.tif")
		if not os.path.exists(ele_tif):
			raise FileNotFoundError(f"[Error] Could not find elevation database")
		values = extract_raster_values(ele_tif, coords)
		df["elevation"] = values
		print(f"[Info] Elevation data extracted")

	elif env_type.lower() == "tavg":
		for month in range(1, 13):
			tavg_tif = os.path.join(data_dir, f"wc2.1_30s_tavg_{month:02d}.tif")
			if not os.path.exists(tavg_tif):
				print(f"[Warning] Could not find average temperature {month} data")
				continue
			values = extract_raster_values(tavg_tif, coords)
			df[f"tavg_{month}"] = values
			print(f"[Info] Average temperature for month{month} extracted")

	elif env_type.lower() == "wind":
		for month in range(1, 13):
			wind_tif = os.path.join(data_dir, f"wc2.1_30s_wind_{month:02d}.tif")
			if not os.path.exists(wind_tif):
				print(f"[Warning] Could not find wind speed {month} data")
				continue
			values = extract_raster_values(wind_tif, coords)
			df[f"wind_{month}"] = values
			print(f"[Info] Wind speed for month{month} extracted")

	elif env_type.lower() == "vapr":
		for month in range(1, 13):
			vapr_tif = os.path.join(data_dir, f"wc2.1_30s_vapr_{month:02d}.tif")
			if not os.path.exists(vapr_tif):
				print(f"[Warning] Could not find water vapor pressure {month} data")
				continue
			values = extract_raster_values(vapr_tif, coords)
			df[f"vapr_{month}"] = values
			print(f"[Info] Water vapor pressure for month{month} extracted")

	elif env_type.lower() == "prec":
		for month in range(1, 13):
			prec_tif = os.path.join(data_dir, f"wc2.1_30s_prec_{month:02d}.tif")
			if not os.path.exists(prec_tif):
				print(f"[Warning] Could not find precipitation {month} data")
				continue
			values = extract_raster_values(prec_tif, coords)
			df[f"prec_{month}"] = values
			print(f"[Info] Precipitation for month{month} extracted")

	elif env_type.lower() == "srad":
		for month in range(1, 13):
			srad_tif = os.path.join(data_dir, f"wc2.1_30s_srad_{month:02d}.tif")
			if not os.path.exists(srad_tif):
				print(f"[Warning] Could not find solar radiation {month} data")
				continue
			values = extract_raster_values(srad_tif, coords)
			df[f"srad_{month}"] = values
			print(f"[Info] Solar radiation for month{month} extracted")

	return df

def extract_all_envs(lon_lat_data, database_dir):
	all_env_data = pd.read_csv(lon_lat_data, sep=',')
	all_env_data.columns = (all_env_data.columns.astype(str)
							.str.replace("\ufeff", "", regex=False)
							.str.strip()
							.str.lower())
	for env_type in ["srad"]:
		print(f"[Info] Extracting {env_type} data...")
		env_data = extracting_envs(lon_lat_data, env_type, database_dir)
		all_env_data = pd.concat([all_env_data, env_data.iloc[:, 2:]], axis=1)

	return all_env_data

if __name__ == "__main__":
	lon_lat_data = "/Users/sy/PhyloSelect/PhyloSelect-v2.1.1/data_creating/geo.csv"
	databse_dir = "/Volumes/T7 Shield/envs"
	result_df = extract_all_envs(lon_lat_data, databse_dir)
	result_df.to_csv("/Users/sy/PhyloSelect/PhyloSelect-v2.1.1/data_creating/geo_envs.csv", index=False)



