import requests
from bs4 import BeautifulSoup
import pandas as pd
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Define the URL of the form
url = 'https://www.camp.bicnirrh.res.in/predict/hii.php'

# Define the file path of the text file you want to upload
file_path = 'C:\\Users\\legion\\Desktop\\evolucijsko\\testni.txt'

# Define the data to be sent with the form
data = {
    'dataset': 'synthetic',  # Change to 'natural' if needed
    'algo[]': ['rf'],  # Select Random Forest
}

# Send POST request to submit the form
with open(file_path, 'rb') as file:
    files = {'userfile': (file.name, file, 'text/plain')}
    response = requests.post(url, data=data, files=files, verify=False)

# Check if the request was successful
if response.status_code == 200:
    # Parse the HTML response
    soup = BeautifulSoup(response.text, 'html.parser')

    # Find the table in the HTML
    table = soup.find('table', attrs={'border': '1'})

    # Get the table rows
    rows = table.find_all('tr')

    # Prepare a list to store the row data
    data = []

    # Iterate over the rows
    for row in rows[1:]:  # Skip the header row
        # Get the columns in each row
        cols = row.find_all('td')
        # Get the text from the columns and add it to the list
        data.append([col.text for col in cols])

    # Create a DataFrame from the list
    df = pd.DataFrame(data, columns=['Seq. ID.', 'Class', 'AMP Probability'])

    # Save the DataFrame to an Excel file
    df.to_excel('C:\\Users\\legion\\Desktop\\output.xlsx', index=False)
else:
    print("Failed to submit the form. Status code:", response.status_code)