import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

st.title("SARS-CoV-2 Mutation Frequency Explorer")

st.markdown("""Upload the processed CSV produced by `build_trends.py` (`mutation_freq_by_day.csv`).""")
uploaded = st.file_uploader("Upload mutation_freq_by_day.csv", type=["csv"])

if uploaded:
    df = pd.read_csv(uploaded, parse_dates=["date"])
    labels = sorted(df["label"].unique())
    selection = st.multiselect("Select mutations to plot", labels[:10], default=labels[:3] if len(labels)>=3 else labels)
    if selection:
        plt.figure()
        for lab in selection:
            sub = df[df["label"] == lab]
            plt.plot(sub["date"], sub["freq"], label=lab)
        plt.xlabel("Date")
        plt.ylabel("Frequency")
        plt.title("Mutation frequencies over time")
        plt.legend()
        st.pyplot(plt.gcf())
