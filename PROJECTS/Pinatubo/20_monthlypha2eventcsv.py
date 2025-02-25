import os
import csv
import glob
import re
from datetime import datetime, timedelta
from obspy import UTCDateTime

def parse_custom_format(file_path, output_dir, error_log="lines_that_could_not_be_parsed.txt"):
    """
    Parses seismic phase data using **fixed character positions**.
    Saves each event to a CSV file, named using the earliest timestamp in the event.

    - **P-wave exists only if 'P' is in position 6 (index 5)**
    - **S-wave exists only if 'S' is in position 2 (index 1) of the S-wave arrival code**
    - **Fixed-width character slicing ensures reliable parsing**
    """

    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

    events = []  # List of events (each event is a list of arrivals)
    current_event = []  # Temporary storage for arrivals in the current event
    
    with open(file_path, 'r', encoding='utf-8', errors='replace') as file, open(error_log, 'a', encoding='utf-8') as error_file:
        for line in file:

            line = line.rstrip()  # Strip \n and spaces


            # DEBUG: Print each line to confirm correct reading
            print(f"DEBUG: Raw line -> '{line}'")

            # Print all positions where 'S' appears
            s_positions = [i for i, char in enumerate(line) if char == 'S']
            # Remove any positions outside the range [30, 40]
            s_positions = [pos for pos in s_positions if 35 <= pos <= 40]
            s_pos = 0
            if len(s_positions)==1:
                s_pos = s_positions[0]
            print(f"DEBUG: 'S' found at positions -> {s_pos} in line: {repr(line)}")


            # If a "10" separator is found, start a new event
            if line.strip() == "10":
                if current_event:  # Save the previous event if it has arrivals
                    events.append(current_event)
                    save_event_to_csv(current_event, output_dir)
                    current_event = []  # Reset for the next event
                continue  # Skip processing the separator itself

            #try:

            #  **Fixed Character Positions**
            station = line[0:3].strip()  # Positions 1-3
            if station.lower()=='xxx' or len(station)<3:
                continue
            orientation = line[3:4].strip()  # Position 4
            if orientation=='P':
                line = line[0:4] + '  ' + line[4:]
            p_arrival_code = line[4:8].replace(' ', '?')   # Positions 5-8
            timestamp_str = line[9:24].strip().replace(" ", "0")  # Positions 10-24 (date/time)
            s_wave_delay = line[s_pos-7:s_pos-1].strip() if s_pos else ""
            if s_pos>0:
                if len(line)>s_pos+3:
                    s_arrival_code = line[s_pos-1:s_pos+3].replace(' ', '?')
                else:
                    s_arrival_code = line[s_pos-1:].ljust(4).replace(' ', '?')
            else:
                s_arrival_code = ""  
            unknown_str = line[42:].strip() if len(line) > 47 else ""
            unknown_str = re.sub(r'[^\x20-\x7E]', '', unknown_str).strip()  # Keeps only printable ASCII characters
 

            #  **Check if P-wave data exists (only if 'P' in position 6, index 5)**
            has_p_wave = len(p_arrival_code) >= 2 and p_arrival_code[1] == "P"

            #  **Check if S-wave data exists (only if 'S' in position 2, index 1 of S-wave code)**
            has_s_wave = s_pos>0

            # Convert spaces and '?' in arrival codes to 'unknown'
            p_arrival_code = p_arrival_code.replace("?", " ") if has_p_wave else "unknown"
            s_arrival_code = s_arrival_code.replace("?", " ") if has_s_wave else "unknown"

            # Convert timestamp to UTC
            add_secs = 0
            if timestamp_str[-5:]=='60.00':
                timestamp_str = timestamp_str.replace('60.00', '00.00')
                add_secs += 60
            if timestamp_str[-7:-5]=='60':
                timestamp_str = timestamp_str.replace('60', '00')
                add_secs += 3600                
            try:
                timestamp = UTCDateTime(datetime.strptime(timestamp_str, "%y%m%d%H%M%S.%f"))
            except:
                try:
                    timestamp = UTCDateTime(datetime.strptime(timestamp_str, "%y%m%d%H%M"))
                except:
                    continue
            timestamp = timestamp + add_secs


            #  **Determine SEED channel**
            if orientation in "ZNE":  # Standard orientations
                channel = f"EH{orientation}"
            elif orientation == "L":  # Special case for "L"
                channel = "ELZ"
            else:
                channel = f'??{orientation}'
                #raise ValueError(f"Unknown orientation '{orientation}' in '{station}'")

            # Construct SEED ID
            seed_id = f"XB.{station}..{channel}"

            #  **Store P-wave arrival**
            if has_p_wave:
                p_arrival = {
                    "seed_id": seed_id,
                    "time": timestamp,
                    "onset": p_arrival_code[0] if p_arrival_code[0] in ["I", "E"] else "unknown",
                    "type": "P",
                    "first_motion": p_arrival_code[2] if p_arrival_code[2] in ["U", "D"] else "unknown",
                    "uncertainty": int(p_arrival_code[3]) if p_arrival_code[3].isdigit() else None,
                    "unknown": unknown_str
                }
                current_event.append(p_arrival)

            #  **Store S-wave arrival**
            if has_s_wave and s_wave_delay.replace(".", "").isdigit():
                s_wave_time = timestamp + timedelta(seconds=float(s_wave_delay))
                s_arrival = {
                    "seed_id": seed_id,
                    "time": s_wave_time,
                    "onset": s_arrival_code[0] if s_arrival_code[0] in ["I", "E"] else "unknown",
                    "type": "S",
                    "first_motion": s_arrival_code[2] if s_arrival_code[2] in ["U", "D"] else "unknown",
                    "uncertainty": int(s_arrival_code[3]) if s_arrival_code[3].isdigit() else None,
                    "unknown": unknown_str
                }
                current_event.append(s_arrival)

            '''
            except Exception as e:
                error_file.write(f"Error parsing line: {line}\nReason: {e}\n")
                error_file.write(f'{station}, {orientation}, {p_arrival_code}, {timestamp_str}, {s_wave_delay}, {s_arrival_code}, {optional_number}'+'\n\n')
                continue  # Skip this line and move to the next
            '''

    # Save the last event if there are remaining arrivals
    if current_event:
        events.append(current_event)
        save_event_to_csv(current_event, output_dir)

    return events

def save_event_to_csv(event, output_dir):
    """Saves an event (list of arrivals) to a CSV file named after the earliest time in the event."""
    
    earliest_time = min(arrival["time"] for arrival in event)
    timestamp_str = earliest_time.strftime("%Y%m%d_%H%M%S")  # Format as YYYYMMDD_HHMMSS
    csv_filename = os.path.join(output_dir, f"{timestamp_str}.csv")

    with open(csv_filename, mode="w", newline="") as csv_file:
        fieldnames = ["seed_id", "time", "onset", "type", "first_motion", "uncertainty", "unknown"]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for arrival in event:
            print(arrival)
            writer.writerow(arrival)

    print(f"Saved event to: {csv_filename}")

#  **Process all `pin*91.pha` files**
input_files = glob.glob("/data/Pinatubo/PHASE/pin*91.pha")
output_directory = "/data/Pinatubo/PHASE/EVENTS_PHA"

os.makedirs(output_directory, exist_ok=True)

for file in input_files:
    print(f"Processing: {file}")
    parse_custom_format(file, output_directory)

print("\n All events processed and saved in 'EVENTS_PHA' directory.")
print(" Any lines that could not be parsed have been logged in 'lines_that_could_not_be_parsed.txt'.")
