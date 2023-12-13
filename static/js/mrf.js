document.getElementById("submit-btn").addEventListener("click", function(event) {
    event.preventDefault();

    var formData = new FormData();
    formData.append("file", document.getElementById("file").files[0]);

    axios.post("/mrf", formData)
        .then(function(response) {
            console.log("Response data:", response.data);
            // Show the results section
            document.getElementById("results").style.display = "block";
            // Show the download link and set its href attribute
            document.getElementById("download-link").style.display = "block";
            document.getElementById("download-link").href = response.data.data;
        })
        .catch(function(error) {
            console.error("Error:", error);
        });
});
