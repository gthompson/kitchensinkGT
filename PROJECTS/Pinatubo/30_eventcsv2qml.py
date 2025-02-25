import os
import csv
import glob
from obspy import UTCDateTime
from obspy.core.event import Event, Origin, Arrival, Pick, WaveformStreamID, Catalog

def csv_to_obspy_event(csv_file):
    """
    Converts an event CSV file into an ObsPy Event object.

    Parameters:
        csv_file (str): Path to the CSV file.

    Returns:
        obspy.core.event.Event: The ObsPy Event object.
    """
    event = Event()  # Create an empty ObsPy Event

    picks = []  # List to store Pick objects
    arrivals = []  # List to store Arrival objects

    with open(csv_file, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file)

        for row in reader:
            # Extract values from CSV
            seed_id = row['seed_id']
            pick_time = UTCDateTime(row['time'])
            onset = row['onset'] if row['onset'] != "unknown" else None
            if onset=='E':
                onset = 'emergent'
            elif onset=='I':
                onset='impulsive'
            else:
                onset=='questionable'
            phase_hint = row['type']  # "P" or "S"
            polarity = row['first_motion'] if row['first_motion'] != "unknown" else None
            if polarity=='U':
                polarity='positive'
            elif polarity=='D':
                polarity='negative'
            else:
                polarity='undecidable'
            uncertainty = float(row['uncertainty']) if row['uncertainty'] else None

            # Create a Pick object
            pick = Pick(
                time=pick_time,
                waveform_id=WaveformStreamID(seed_string=seed_id),
                onset=onset,
                phase_hint=phase_hint,
                polarity=polarity
            )
            picks.append(pick)

            # Create an Arrival object (linking to the Pick)
            arrival = Arrival(
                pick_id=pick.resource_id,
                phase=phase_hint,
                time_residual=uncertainty if uncertainty is not None else 0.0
            )
            arrivals.append(arrival)

    # Create an Origin object (hypocenter)
    origin_time = min(pick.time for pick in picks)  # Use earliest pick as origin time
    origin = Origin(
        time=origin_time,
        arrivals=arrivals
    )

    # Add Origin and Picks to Event
    event.origins.append(origin)
    event.picks.extend(picks)

    return event

def save_catalog_as_quakeml(catalog, output_file):
    """
    Saves an ObsPy Catalog object as a single QuakeML file.

    Parameters:
        catalog (obspy.core.event.Catalog): The ObsPy Catalog containing multiple events.
        output_file (str): Path to save the QuakeML file.
    """
    catalog.write(output_file, format="QUAKEML")
    print(f"✅ QuakeML file saved: {output_file}")

# 🔹 **Process all event CSV files**
output_directory = "/data/Pinatubo/PHASE/EVENTS_PHA"
csv_files = glob.glob("/data/Pinatubo/PHASE/EVENTS_PHA/*.csv")  # Adjust the folder path if needed
quakeml_output = "/data/Pinatubo/PHASE/all_events.xml"  # Name of the output QuakeML file

# Create an empty ObsPy Catalog
catalog = Catalog()

for csv_file in csv_files:
    print(f"Processing: {csv_file}")
    event = csv_to_obspy_event(csv_file)  # Convert CSV to Event
    catalog.events.append(event)  # Add to Catalog

# Save the full Catalog as a QuakeML file
save_catalog_as_quakeml(catalog, quakeml_output)

print("\n✅ All events converted and saved in 'all_events.xml'!")
