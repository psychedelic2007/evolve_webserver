document.getElementById('uploadForm').addEventListener('submit', function(event) {
    event.preventDefault(); // Prevent the default form submission

    // Get form data
    var formData = new FormData(event.target);

    // Send the form data to the server using Fetch API
    fetch('/mi_calculate', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        if (data.error_message) {
            alert(data.error_message); // Display an alert for any error message
        } else if (data.graph_data) {
            // Update the image source
            document.getElementById('networkPlot').src = 'data:image/png;base64,' + data.graph_data;

            // Display the download link
            var downloadLink = document.getElementById('downloadLink');
            downloadLink.href = data.processed_data_link;
            downloadLink.style.display = 'block';
        }
    })
    .catch(error => {
        console.error('Error:', error);
        alert('An error occurred while processing the file.'); // Display an error alert
    });
});
