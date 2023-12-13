@app.route('/prediction', methods=['GET','POST'])
def preediction_try():
    if session.get('logged_in'):
        user = User.query.filter_by(username=session['username']).first()
        return render_template('prediction.html', user=user, num_files=1)
    else:
        return redirect(url_for('login'))

@app.route('/predict', methods=['POST'])
def prediction():
    if request.method == 'POST':
        file = request.files['file']

        if file:
            data = pd.read_csv(file)
            x = data.drop('Target', axis=1)
            # Load your model (replace with your model loading code)
            model = load_model("model_final.h5")
            prediction = model.predict(x)

            # Save prediction data to a CSV file
            prediction_file = pd.DataFrame({
                "Residue Position": range(1, len(data) + 1),
                "Mutational Probability": prediction.flatten()
            })
            prediction_csv_path = "static/prediction_data.csv"
            prediction_file.to_csv(prediction_csv_path, index=False)

            # Generate the plot
            plt.plot(prediction)
            plt.xlabel("Residue Position")
            plt.ylabel("Mutational Entropy")
            plt.title("Prediction Plot")
            img_buffer = io.BytesIO()
            plt.savefig(img_buffer, format="png")
            img_buffer.seek(0)
            img_data = base64.b64encode(img_buffer.read()).decode()
            plt.close()

        if session.get('logged_in'):
            user = User.query.filter_by(username=session['username']).first()
            return render_template('prediction.html', user=user, graph_data=img_data)
        else:
            return redirect(url_for('login'))
    else:
        return render_template('index.html')

@app.route('/download_prediction_csv')
def download_prediction_csv():
    # Provide a link to download the prediction data as a CSV file
    prediction_csv_path = "static/prediction_data.csv"
    return send_file(prediction_csv_path, as_attachment=True, download_name="prediction_data.csv")

@app.route('/download_prediction_plot')
def download_prediction_plot():
    plot_path = "static/prediction_plot.png"
    import os
    if os.path.exists(plot_path):
        try:
            return send_file(plot_path, as_attachment=True, download_name="prediction_plot.png")
        except Exception as e:
            print(f"Error while sending the image: {str(e)}")
            return "An error occurred while sending the image."
    else:
        return "The prediction plot is not available."