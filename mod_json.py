import json
import os


def modify_json(json_file):
    """
    Modify the JSON file to rescale the TLS annotations. Saves the modified
    JSON file to the mod_TLS_json directory as the original file and returns the modified
    JSON object.

    Args:
        json_file (dict): The JSON file to modify

    Returns:
        dict: The modified JSON file
    """
    # Open the JSON file
    with open(json_file, 'r') as file:
        json_data = json.load(file)

    # Rescale the TLS annotations
    for feature in json_data['features']:
        for polygon in feature['geometry']['coordinates']:
            for coord in polygon:
                for i in range(len(coord)):
                    coord[i] *= 2
    
    base_name = os.path.basename(json_file)
    mod_file_name = f"mod_{base_name}"
    output_path = os.path.join('mod_TLS_json', mod_file_name)
    
    _save_json(json_data, output_path)

    return json_data


def _save_json(json_data, output_file):
    """
    Save the modified JSON data to a new file.

    Args:
        json_data (dict): The modified JSON data
        output_file (str): The path to the output file
    """

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Save the modified JSON data to the output file
    with open(output_file, 'w') as outfile:
        json.dump(json_data, outfile, indent=4)